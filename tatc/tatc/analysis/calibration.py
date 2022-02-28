# add cross-calibration analysis functions
import orekit
import math
from datetime import datetime, timedelta
vm = orekit.initVM()
import pandas as pd
from orekit.pyhelpers import setup_orekit_curdir
from pkg_resources import resource_stream
setup_orekit_curdir()
from tatc_core.generation.points import generate_fibonacci_lattice_points
from tatc_core.generation.cells import generate_cubed_sphere_cells
from org.orekit.orbits import KeplerianOrbit, PositionAngle
from org.orekit.bodies import CelestialBodyFactory
from org.orekit.frames import FramesFactory
from org.orekit.utils import Constants, PVCoordinatesProvider
from org.orekit.propagation.analytical.tle import TLE
from orekit.pyhelpers import datetime_to_absolutedate
import pandas as pd
import geopandas as gpd
import geoplot as gplt
import numpy as np

import time
import math
from datetime import datetime, timedelta, timezone
from joblib import Parallel, delayed
from shapely.geometry import box, Polygon

def continent_filter(obsevation_sample_distance,po):

    #obsevation_sample_distance = 5e5

    df=pd.DataFrame(columns=["id","geometry"])
    filter_po=[]
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    for j in range(len(po)):
        for i in range(len(world["geometry"])):
            if po["geometry"][j].within(world["geometry"][i])==True:
                filter_po.append(po["id"][j]) #gets point ids that are located in the continent

    for i in range(len(po)):
        re= any(po["id"][i]== ids for ids in filter_po)
        if re==False:
            df=df.append({"id":po["id"][i],"geometry":po["geometry"][i]},ignore_index=True)
    points=  gpd.GeoDataFrame(
          df)
    base = world.plot(color='white', edgecolor='black')
    points.plot(ax=base, marker='o', color='red', markersize=1)
    return points

def generate_walker_delta_members(constellation):
    """
    Generate member satellites of a walker-delta constellation

    Args:
        constellation (:obj:`WalkerConstellation`): The walker-delta constellation object.
    Returns:
        satellites (:obj:List) The list of member satellites.
    """

    satellites = []
    for satellite in range(constellation["number_satellites"]):
        satellites_per_plane = math.ceil(constellation["number_satellites"]/constellation["number_planes"])
        plane = satellite // satellites_per_plane
        num_planes = constellation["number_planes"]
        rel_spacing = constellation["relative_spacing"]
        lead_orbit = constellation["lead_orbit"]
        lead_tle = TLE(lead_orbit["tle"][0], lead_orbit["tle"][1])
        tle = TLE(
            0, # satellite number
            'U', # classification
            0, # launch year
            0, # launch number
            'A  ', # launch piece
            0, # ephemeris type
            0, # element number
            lead_tle.getDate(), # epoch
            lead_tle.getMeanMotion(), # mean motion
            0.0, # mean motion first derivative
            0.0, # mean motion second derivative
            lead_tle.getE(), # eccentricity
            lead_tle.getI(), # inclination
            lead_tle.getPerigeeArgument(), # periapsis argument
            lead_tle.getRaan() + plane*2*math.pi/num_planes, # right ascension of ascending node
            # FIXME orekit constructor requires mean anomaly, rather than true anomaly (as specified for walker delta constellations)
            lead_tle.getMeanAnomaly() + ((satellite % satellites_per_plane)*num_planes
                         + rel_spacing*plane)*2*math.pi/(satellites_per_plane*num_planes), # mean anomaly
            0, # revolution number at epoch
            0.0 # b-star (ballistic coefficient)
        )
        satellites.append(
            {"name" : constellation["name"] + "{:03d}".format(satellite+1),
            "orbit" : [
                tle.getLine1(),
                tle.getLine2()
            ]}

        )

    return satellites

class obs():
    """
    This class determines the start and end time observations
    obs: observations
    """
    def __init__(self, name ,start, end, sun):
        self.name=name
        self.start= start
        self.end=end
        self.sun=sun

def cross_opportunities(df,point_ids,max_range):
    overlaps=[]
    i=0

    for j in range(len(point_ids)):

        while point_ids[j]==df['id'][i]:
            if df['satellite'][i]!= "Aqua":
                if df['satellite'][i+1]=="Aqua":
                    r1= obs(df['satellite'][i],df['start'][i],end=df['end'][i], sun=df["sat_sunlit"][i])
                    r2= obs(df['satellite'][i+1],df['start'][i+1],end=df['end'][i+1], sun=df["sat_sunlit"][i+1])
                    if r2.start<= r1.end and r2.end>=r1.end:
                        delta= (r1.end-r2.start)/timedelta(seconds=1)
                        overlaps.append({"overlap_time": delta,
                                         "type": "simultaneous",
                                         "name": r1.name,
                                         "leading_sat": r1.name,
                                         "trailing_sat": r2.name,
                                         "start_time_leading": r1.start,
                                         "start_time_trailing": r2.start,
                                         "end_time_leading": r1.end,
                                         "end_time_trailing": r2.end,
                                         "start_obs": r2.start,
                                         "end_obs": r1.end,
                                         "sun_light":r1.sun,
                                         'point_id': point_ids[j]})

                    elif r2.start<= r1.end and r2.end<r1.end:
                        delta= (r2.end-r2.start)/timedelta(seconds=1)
                        overlaps.append({"overlap_time": delta,
                                         "type": "simultaneous",
                                         "name": r1.name,
                                         "leading_sat": r1.name,
                                         "trailing_sat": r2.name,
                                         "start_time_leading": r1.start,
                                         "start_time_trailing": r2.start,
                                         "end_time_trailing": r2.end,
                                         "end_time_leading": r1.end,
                                         "start_obs": r2.start,
                                         "end_obs": r2.end,
                                         "sun_light":r1.sun,
                                         'point_id': point_ids[j]})

                    elif r2.start>=r1.end and r2.end<= r1.end+max_range:
                        delta=(r2.end-r2.start)/timedelta(seconds=1)
                        overlaps.append({"overlap_time": delta,
                                         "type": "pseudo_invariant",
                                         "name": r1.name,
                                         "leading_sat": r1.name,
                                         "trailing_sat": r2.name,
                                         "start_time_leading": r1.start,
                                         "start_time_trailing": r2.start,
                                         "end_time_leading": r1.end,
                                         "end_time_trailing": r2.end,
                                         "end_time_window": r1.end+max_range,
                                         "start_obs":r2.start,
                                         "end_obs": r2.end,
                                         "sun_light":r1.sun,
                                         'point_id': point_ids[j]})

                    elif r2.start<r1.end+max_range and r2.end>=r1.end+max_range:
                        delta= (r1.end+ max_range-r2.start)/timedelta(seconds=1)
                        overlaps.append({"overlap_time": delta,
                                         "type": "pseudo_invariant",
                                         "name": r1.name,
                                         "leading_sat": r1.name,
                                         "trailing_sat": r2.name,
                                         "start_time_leading": r1.start,
                                         "start_time_trailing": r2.start,
                                         "end_time_leading": r1.end,
                                         "end_time_trailing": r2.end,
                                         "end_time_window": r1.end+max_range,
                                         "start_obs": r2.start,
                                         "end_obs": r1.end+max_range,
                                         "sun_light":r1.sun,
                                         'point_id': point_ids[j]})


            elif df['satellite'][i]== "Aqua":
                if df['satellite'][i+1]!= "Aqua":
                    r1= obs(df['satellite'][i],df['start'][i],end=df['end'][i],sun=df["sat_sunlit"][i])
                    r2= obs(df['satellite'][i+1],df['start'][i+1],end=df['end'][i+1], sun=df["sat_sunlit"][i+1])
                    if r2.start<= r1.end and r2.end>=r1.end:
                        delta= (r1.end-r2.start)/timedelta(seconds=1)
                        overlaps.append({"overlap_time": delta,
                                         "type":"simultaneous",
                                         "name": r2.name,
                                         "leading_sat": r1.name,
                                         "trailing_sat": r2.name,
                                         "start_time_leading": r1.start,
                                         "start_time_trailing": r2.start,
                                         "end_time_leading": r1.end,
                                         "end_time_trailing": r2.end,
                                         "start_obs": r2.start,
                                         "end_obs": r1.end,
                                         "sun_light": r2.sun,
                                         'point_id': point_ids[j]})

                    elif r2.start<= r1.end and r2.end<=r1.end:
                        delta= (r2.end-r2.start)/timedelta(seconds=1)
                        overlaps.append({"overlap_time": delta,
                                         "type":"simultaneous",
                                         "name": r2.name,
                                         "leading_sat": r1.name,
                                         "trailing_sat": r2.name,
                                         "start_time_leading": r1.start,
                                         "start_time_trailing": r2.start,
                                         "end_time_trailing": r2.end,
                                         "end_time_leading": r1.end,
                                         "start_obs": r2.start,
                                         "end_obs": r2.end,
                                         "sun_light": r2.sun,
                                         'point_id': point_ids[j]})

                    elif r2.start>=r1.end and r2.end<= r1.end+max_range:
                        delta=(r2.end-r2.start)/timedelta(seconds=1)
                        overlaps.append({"overlap_time": delta, #units of sec
                                         "type": "pseudo_invariant",
                                         "name": r2.name,
                                         "leading_sat": r1.name,
                                         "trailing_sat": r2.name,
                                         "start_time_leading": r1.start,
                                         "start_time_trailing": r2.start,
                                         "end_time_leading": r1.end,
                                         "end_time_trailing": r2.end,
                                         "end_time_window": r1.end+max_range,
                                         "start_obs":r2.start,
                                         "end_obs": r2.end,
                                         "sun_light": r2.sun,
                                         'point_id': point_ids[j]})

                    elif r2.start<r1.end+max_range and r2.end>=r1.end+max_range:
                        delta= (r1.end+ max_range-r2.start)/timedelta(seconds=1)
                        overlaps.append({"overlap_time": delta,
                                         "type": "pseudo_invariant",
                                         "name": r2.name,
                                         "leading_sat": r1.name,
                                         "trailing_sat": r2.name,
                                         "start_time_leading": r1.start,
                                         "start_time_trailing": r2.start,
                                         "end_time_leading": r1.end,
                                         "end_time_trailing": r2.end,
                                         "end_time_window": r1.end+max_range,
                                         "start_obs": r2.start,
                                         "end_obs": r1.end+max_range,
                                         "sun_light": r2.sun,
                                         'point_id': point_ids[j]})


            if i< len(df)-2:
                i+=1
            else:
                i=0

    return overlaps

import numpy as np
import pandas as pd
from datetime import timedelta

def average_calibration(df):
    id_list=set(df['id'])
    names_list= set(df['satellite'])
    inst_list= set(df['instrument'])
    sat_names=[]
    point_ids=[]
    instrument_list=[]
    for i in id_list:
        point_ids.append(i)
    for i in names_list:
        sat_names.append(i)
    for i in inst_list:
        instrument_list.append(i)
    point_ids.sort()    # sort ids from minimum to maximum


    max_range=timedelta(hours=1)
    detectors= cross_opportunities(df,point_ids,max_range)

    aver= []
    for i in range(len(detectors)):
        aver.append(detectors[i]["overlap_time"])

    res= np.average(aver)

    data_point=[]
    sat_opp=[]

    for i in range (len(point_ids)):
        obser=[]
        for j in range(len(detectors)):
            if point_ids[i]== detectors[j]["point_id"]:
                obser.append({"start":detectors[j]["start_obs"],
                            "end": detectors[j]["end_obs"]})
        if obser!=[]:
            data_point.append({"id":point_ids[i],
                            "timeline":obser})
    sat_cross=[]
    for sat in sat_names:
        sam=[]
        opp=pd.DataFrame(columns=["start","end"])
        for j in range(len(detectors)):
            if sat== detectors[j]["name"]:
                if detectors[j]["sun_light"]==True:
                    opp=opp.append({"start": detectors[j]["start_obs"], "end": detectors[j]["end_obs"]},ignore_index=True)
        rev_persat=[]
        opp=opp.sort_values("start").reset_index(drop=True)

        for i in range(len(opp)):
            sam.append({
                            "start":opp["start"][i],
                            "end": opp["end"][i]
                            })
        a=0
        max_len= len(sam)
        add=[]
        for i in range(len(sam)):

            if len(add)==0:
                add.append(sam[i])
                max_tim=sam[i]["end"]
            elif sam[i]["start"]<= max_tim:
                a+=1
            elif sam[i]["start"]> max_tim:
                add.append(sam[i])
                max_tim= sam[i]["end"]

        sat_cross.append({"name": sat,
                        "timeline_opp": add})
        del opp["start"]
        del opp["end"]


    #Compute average revisit time per satellite
    all_sat=[]
    sat_av=0
    for i in range(len(sat_cross)):
        rev_list=[]
        for j in range(len(sat_cross[i]["timeline_opp"])-1):
            a=(sat_cross[i]["timeline_opp"][j+1]["start"]- sat_cross[i]["timeline_opp"][j]["end"])/timedelta(hours=1)
            if a> 24: # values greater than 24 hours for cross calibration opportunities
                rev_list.append(a)

        if sat_cross[i]["name"]!= "Aqua":
            all_sat.append({
                "sat_name": sat_cross[i]["name"],
                "revisit_opp": np.average(rev_list)
                })
            sat_av+=np.average(rev_list)

    return  sat_av/(len(sat_names)-1)  ## Units in hours

def cross_calibration(sat_min,sat_max,duration):
    mission_duration = timedelta(days=duration)
    instrument = Instrument(
    'MODIS',
    field_of_regard=115.39, # degrees, determines point observability
    req_access_time=timedelta(seconds=10) # determines valid observations
    )
    inst_Tesat=Instrument(
    'Planet',
    field_of_regard= 15,  #Flock instrument
    req_access_time=timedelta(seconds=10)
    )

    ##### defines grid points  #####
    target_region = gpd.GeoDataFrame({'geometry': [box(-160, 50, -50, 23)]}, crs="EPSG:4326")
    obsevation_sample_distance = 5e5
    grid_sample_distance = 0.5e6
    po = generate_fibonacci_lattice_points(obsevation_sample_distance, target_region)
    points= continent_filter(obsevation_sample_distance,po)
    grid_cells = generate_cubed_sphere_cells(grid_sample_distance, target_region, strips='lat')
    #########################################
    ################## Aqua TLEs  ##################
    tle1='1 27424U 02022A   21069.59193972  .00000087  00000-0  29406-4 0  9998'
    tle2='2 27424  98.2123  12.3766 0000444 141.7042  31.1842 14.57113601  2631'
    ################################
    fin_result= pd.DataFrame(columns=["#sat", "calibration_opp_hr","revisit_hr"])
    ####################ISS Constellation######
    for n_sat in range(sat_min, sat_max):
        constellation= {
                "name": "Test-one",
                "lead_orbit": {
                "type": "tle",
                "tle": [
                    "1 25544U 98067A   21133.52830052  .00001043  00000-0  27100-4 0  9997",
                    "2 25544  51.6441 154.1527 0003110 356.5241 135.9982 15.48997083283227"
                ]
                },
                "number_satellites": n_sat,
                "number_planes": 4,
                "relative_spacing": 15
            }
        a= generate_walker_delta_members(constellation)

        df=pd.DataFrame(columns=["name","tle1","tle2"])
        for i in range(len(a)):
            df= df.append({"name":a[i]["name"], "tle1": a[i]["orbit"][0], "tle2": a[i]["orbit"][1]},ignore_index=True)
        satellites=[
            Satellite(df['name'][i], tle=[df["tle1"][i], df["tle2"][i]],instruments=[inst_Tesat])
            for i in range(len(df))
            ]
        satellites.append(Satellite('Aqua', tle=[tle1, tle2],instruments=[instrument]))

        # collect raw observations

        tic = time.perf_counter()

        raw_observations = pd.concat(
            Parallel(n_jobs=-1)(
                delayed(collect_observations)(
                    Point(r.id, r.geometry.y, r.geometry.x),
                    satellite,
                    instrument,
                    datetime.now(timezone.utc),
                    datetime.now(timezone.utc) + mission_duration
                )
                for satellite in satellites
                for instrument in satellite.instruments
                for _, r in points.iterrows()
            ),
            ignore_index=True
        )
        toc = time.perf_counter()
        print(f"Collection completed in {toc-tic:0.4f} seconds")
        # reduce raw observations to group overlapping observations
        tic = time.perf_counter()

        reduced_observations = pd.concat(
            Parallel(n_jobs=-1)(
                delayed(reduce_observations)(raw_observations.loc[raw_observations['id']==id])
                for id in raw_observations['id'].unique()
            ),
            ignore_index=True
        )

        toc = time.perf_counter()

        print(f"Reduction completed in {toc-tic:0.4f} seconds")
        # plot the access statistic as points
        reduced_point_data, ax = plot_points_access(points, reduced_observations, map_projection)#, extent=(-130,40,-70, 50 ))
        # plot the revisit statistic as points
        reduced_point_data, ax = plot_points_revisit(points, reduced_observations, map_projection)#, extent=(-130,40,-70, 50 ))
        # plot the revisit statistic as grid cells
        grid_data, ax = plot_cells_access(grid_cells, reduced_point_data, gplt.PlateCarree())
        # plot the revisit statistic as grid cells
        grid_data, ax = plot_cells_revisit(grid_cells, reduced_point_data, gplt.PlateCarree())
        # plot observation angles for an example point
        fig = plot_observation_angles(raw_observations.loc[raw_observations['id']==points.iloc[int(len(points)/2)]['id']])


        ####################### writes results files for tracking purposes
        raw_observations.sort_values(by=['id','start'], inplace=True)
        raw_observations.to_csv('raw_observations.csv')
        reduced_observations.to_csv('reduced_observations.csv')
        ##########################
        df= pd.read_csv('raw_observations.csv')
        rev_ti=pd.read_csv('reduced_observations.csv')

        for i in range(len(df['start'])):  ##This loop converts the strings read in the csv file to timestamp format
            df['epoch'][i]=pd._libs.tslibs.timestamps.Timestamp(df['epoch'][i])
            df['end'][i]=pd._libs.tslibs.timestamps.Timestamp(df['end'][i])
            df['start'][i]=pd._libs.tslibs.timestamps.Timestamp(df['start'][i])
        cross_ave= average_calibration(df)


        emp=[]
        for i in range(len(rev_ti['revisit_hr'])):
            if math.isnan(rev_ti['revisit_hr'][i])==False:
                emp.append(rev_ti['revisit_hr'][i])

        revisit_time=np.average(emp)

        fin_result= fin_result.append({"#sat":n_sat,
                                    "calibration_opp_hr":cross_ave,
                                    "revisit_hr":revisit_time},ignore_index=True)
    return fin_results
