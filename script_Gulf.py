#!/usr/bin/env python
"""
"""
from datetime import datetime, timedelta
TOTAL_start_time = datetime.now()

import os
from datetime import datetime, timedelta

from gnome import basic_types
from gnome.spill.elements import ElementType


from gnome import scripting
from gnome.basic_types import datetime_value_2d
from gnome import utilities

from gnome.utilities.remote_data import get_datafile

from gnome.model import Model

from gnome.map import MapFromBNA
from gnome.spill import point_line_release_spill

from gnome.movers import RandomMover, constant_wind_mover, GridCurrentMover
from gnome.movers import GridWindMover

from gnome.outputters import (Renderer,
                               NetCDFOutput
                              )
from gnome.outputters import WeatheringOutput
from gnome.environment import Water, Waves, constant_wind, Wind
from gnome.weatherers import (Emulsification,
                              Evaporation,
                              NaturalDispersion,
                              ChemicalDispersion,
                              Burn,
                              Skimmer,
                              WeatheringData)
                              
from gnome.basic_types import numerical_methods
from gnome.persist import load
from numpy.random import uniform, triangular

import array
import numpy
from netCDF4 import Dataset
from shapely.geometry import Point, Polygon

import matplotlib.pyplot as plt
import time
from pymsgbox import *

# Define base directories------------------------------------------------------
base_dir = os.path.dirname(__file__)
input_dir=os.path.join(base_dir,'input_files')
if not os.path.exists(input_dir):
    os.makedirs(input_dir)
    input_alert = open(input_dir+'/README.txt','w')
    input_alert.write("Four text files are required. You should make these files like the template.")
    input_alert.write("\n1-basic_parameters")
    input_alert.write("\n2-simulation_parameters")
    input_alert.write("\n3-spill_locations")
    input_alert.write("\n4-polygons")

    input_alert.write("\n--------------------------------------------------------------------------")
    input_alert.close()
    alert(text='You have not created the required input text files. Now a folder with name input_files was created in the directory of your python code. Please read the README.txt file in this folder for help.', title='Input Text Files are Required', button='OK')

#----------------------------------------------------------------------------

#load the number of spills and number of start times and number of iterations 
with open(input_dir+'/basic_parameters.txt') as fp:
    for i in fp.readlines():
        data=i.split("-")
        number_of_spills=int(data[1])
        number_of_starttimes=int(data[2])
        number_iterations_for_starttime=int(data[3])
fp.close()
#number of total simulations to be done:
number_of_simulations=number_of_spills*number_of_starttimes*number_iterations_for_starttime
#------------------------------------------------------------------------------

# Load the simulation duration and time between NC snapshots(averaging time)
with open(input_dir+'/simulation_parameters.txt') as fp:
    for i in fp.readlines():
        data=i.split("-")
        days_of_simulation=int(data[1])
        hours_between_snapshots=int(data[2])
fp.close()
simulation_duration=timedelta(days=days_of_simulation)
NC_timestep=timedelta(hours=hours_between_snapshots)
simulation_duration_float=simulation_duration.total_seconds()
NC_timestep_float=NC_timestep.total_seconds()
NC_steps=int(simulation_duration_float/NC_timestep_float)

#------------------------------------------------------------------------------



#Loading Spill locations-------------------------------------------------------
spill_location=[]
oil_type=[]
source_strength=[] # defined to account for the number of wells in each source
with open(input_dir+'/spill_locations.txt') as fp:
    for i in fp.readlines():
        data = i.split("-")
        spill_location.append((float(data[1]), float(data[2]),float(data[3])))
        oil_type.append(str(data[4]))
        source_strength.append(float(data[5]))
fp.close()
#------------------------------------------------------------------------------


# Loading Polygon Zones-------------------------------------------------------
P=[]
sensitivity=[]
with open(input_dir+'/polygons.txt') as fp:
    for i in fp.readlines():
        p1=[]
        data = i.split("-")
        l=len(data)
        for j in range(1,((l-3)/2)+1):
            p1.append((float(data[(j*2)-1]), float(data[j*2])))
        P.append(p1)
        sensitivity.append(float(data[l-2]))
fp.close()
number_of_polygons=len(P)

#------------------------------------------------------------------------------
#Defining matrices spill amount and frequency
start_time_matrix=[[[] for y in xrange(number_of_simulations)] for z in xrange(number_of_spills) ]
gulf_amount_matrix=[[[[0]for x in xrange(number_iterations_for_starttime)]for y in xrange(number_of_starttimes)] for z in xrange(number_of_spills)]
gulf_frequency_matrix=[[[[0]for x in xrange(number_iterations_for_starttime)]for y in xrange(number_of_starttimes)] for z in xrange(number_of_spills)]
gulf_duration_matrix=[[[[0]for x in xrange(number_iterations_for_starttime)]for y in xrange(number_of_starttimes)] for z in xrange(number_of_spills)]

start_time_master_matrix=[[[[0]for x in xrange(number_iterations_for_starttime)]for y in xrange(number_of_starttimes)] for z in xrange(number_of_spills)]

#defining matrices for time random sampling
gulf_month=array.array('i',(0 for i in range (0,number_of_simulations))) # integer only
gulf_day=array.array ('i',(0 for i in range (0,number_of_simulations)))
gulf_hour=array.array ('i',(0 for i in range (0,number_of_simulations)))
#------------------------------------------------------------------------------


dir_gen='Input_generator_result'
if not os.path.exists(dir_gen):
    os.makedirs(dir_gen)
my_input_sets = open(dir_gen+'/README.txt','w')
my_input_sets.write("For each spill location, row number is the start time number and column number is the iteration number")
my_input_sets.write("\n------------------------------------------------------------------------------------------------------")
my_input_sets.close()
# Input generator separated from heart of model
for s in range (0,number_of_spills):
    for i in range(0,number_of_starttimes):
        gulf_month[i]=numpy.floor (uniform(1,13)).astype(int)
        if gulf_month[i]==1 or gulf_month[i]==3 or gulf_month[i]==5 or gulf_month[i]==7 or gulf_month[i]==8 or gulf_month[i]==10:
            gulf_day[i]=numpy.floor (uniform(1,32)).astype(int)
        elif gulf_month[i]==12:
            gulf_day[i]=numpy.floor (uniform (1,(32-(simulation_duration.days+1)))).astype(int) # because of data limitation at the end of time interval
        elif gulf_month[i]==4 or gulf_month[i]==6 or gulf_month[i]==9 or gulf_month[i]==11:
            gulf_day[i]=numpy.floor (uniform(1,31)).astype(int)
        elif gulf_month[i]==2:
            gulf_day[i]=numpy.floor (uniform(1,29)).astype(int)
            
        if not gulf_month[i]==1 and gulf_day[i]==1:
            gulf_hour[i]=numpy.floor(uniform(0,24)).astype(int)
        else:gulf_hour[i]=numpy.floor(uniform(3,24)).astype(int)# my data limitation
        
        start_time= datetime(2015, gulf_month[i], gulf_day[i],gulf_hour[i])
        start_time_matrix[s][i]=start_time
        
        for q in range(0,number_iterations_for_starttime):
            start_time_master_matrix [s][i][q]=start_time_matrix[s][i]
            gulf_amount_matrix[s][i][q]=triangular(50,10000,150000)
            gulf_frequency_matrix[s][i][q]=(6.9*10**-5*(gulf_amount_matrix[s][i][q])**-0.3)*source_strength[s]
            gulf_duration_matrix[s][i][q]=triangular(1,15,30)
        
            

    numpy.savetxt(dir_gen+'/Time_for_Spill%i.txt'%(s),start_time_master_matrix [s][:][:],fmt='%s')
    numpy.savetxt(dir_gen+'/Spill_amount_Spill%i.txt'%(s),gulf_amount_matrix[s][:][:])
    numpy.savetxt(dir_gen+'/Spill_frequency_Spill%i.txt'%(s),gulf_frequency_matrix[s][:][:])
    numpy.savetxt(dir_gen+'/Spill_duration_Spill%i.txt'%(s),gulf_duration_matrix[s][:][:])
        
print '\nThe set of input parameters have been generated and saved successfully'
print '\nNow the simulations will be started'
time.sleep(2)

# running the heart of model------------------------------------------------------
simulation_start_time = datetime.now()

for s in range (0,number_of_spills):
    for i in range(0,number_of_starttimes):
        for q in range(0,number_iterations_for_starttime):
            def make_model():
                
                #Generating water object (on monthly basis)--------------------
                if start_time_master_matrix[s][i][q].month==1: 
                    Temperature=290 
                    Salinity=35
                elif start_time_master_matrix[s][i][q].month==2: 
                    Temperature=295
                    Salinity=35
                elif start_time_master_matrix[s][i][q].month==3: 
                    Temperature=300
                    Salinity=35
                elif start_time_master_matrix[s][i][q].month==4: 
                    Temperature=300
                    Salinity=35
                elif start_time_master_matrix[s][i][q].month==5: 
                    Temperature=305 
                    Salinity=35
                elif start_time_master_matrix[s][i][q].month==6: 
                    Temperature=305
                    Salinity=35
                elif start_time_master_matrix[s][i][q].month==7: 
                    Temperature=310
                    Salinity=35
                elif start_time_master_matrix[s][i][q].month==8: 
                    Temperature=310
                    Salinity=35
                elif start_time_master_matrix[s][i][q].month==9: 
                    Temperature=300
                    Salinity=35
                elif start_time_master_matrix[s][i][q].month==10: 
                    Temperature=300
                    Salinity=35
                elif start_time_master_matrix[s][i][q].month==11: 
                    Temperature=295 
                    Salinity=35
                elif start_time_master_matrix[s][i][q].month==12: 
                    Temperature=290
                    Salinity=35
                    
                water = Water(Temperature, Salinity) #kelvin & PSU
                wind = constant_wind(5, uniform(0,360), 'meters per second') #should be studied
                waves = Waves(wind, water)
                
                # the environment objects created locally (not global)---------
                
                print 'initializing the model'
                model = Model(start_time=start_time_master_matrix[s][i][q], duration=simulation_duration,
                              time_step=4 * 3600, uncertain=True)
                            
                mapfile = get_datafile(os.path.join(base_dir, 'gulf.bna'))
            
                print 'adding the map'
                model.map = MapFromBNA(mapfile, refloat_halflife=6)  # hours
            
                #
                # Add the outputters -- render to images, and save out as netCDF
                #
                root_image_dir=os.path.join(base_dir,'Results','Spill %i'%s,'Images','Start_time %i'%i)
                if not os.path.exists(root_image_dir):
                    os.makedirs(root_image_dir)
                
                print 'adding renderer'
                gulf_images_dir=os.path.join(root_image_dir,'Iteration %i'%q)
                model.outputters += Renderer(mapfile,
                                             gulf_images_dir,
                                             size=(800, 600),)
                                             
                                             #draw_back_to_fore=True)
            
                model.outputters += WeatheringOutput()
            
                
                print "adding netcdf output"
                root_nc_dir=os.path.join(base_dir, 'Results','Spill %i'%s,'NC','Start_time %i'%i)
                if not os.path.exists(root_nc_dir):
                    os.makedirs(root_nc_dir)
                netcdf_output_file = os.path.join(root_nc_dir, 'gulf_output%i.nc'%q)
                scripting.remove_netcdf(netcdf_output_file)
                model.outputters += NetCDFOutput(netcdf_output_file, which_data='all',output_timestep=NC_timestep)
#                
                
                #
                # Set up the movers:
                #
            
                print 'adding a RandomMover:'
                model.movers += RandomMover(diffusion_coef=10000)
            
                print 'adding a simple wind mover:'
            #    model.movers += constant_wind_mover(5, 315, units='m/s')
                model.movers+=GridWindMover('wind_final.nc')
            
                print 'adding a current mover:'
            
        #        # # this is HYCOM currents
                curr_file = get_datafile(os.path.join(base_dir, 'CMEMS_final.nc'))
                model.movers += GridCurrentMover(curr_file,
                                                 num_method=numerical_methods.euler);
            
                # #
                # # Add some spills (sources of elements)
                # #
            
                print 'adding one spill'   
                end_time = start_time_master_matrix[s][i][q] + timedelta(days=gulf_duration_matrix[s][i][q])
                spill = point_line_release_spill(num_elements=1000,
                                                 start_position=spill_location[s],
    #                                             start_position=(55,
    #                                                             26, 0.0),
                                                 release_time=start_time_master_matrix[s][i][q],
                                                 end_release_time=end_time,
    #                                             amount=1000,
                                                 amount=gulf_amount_matrix[s][i][q],
    #                                             substance='NOWRUZ',
                                                 substance=oil_type[s],
    #                                             substance='oil_4',
                                                 units='bbl', water=water)
            
                model.spills+=spill
                
                
#                skim1_start = start_time + timedelta(hours=15.58333)
#                skim2_start = start_time + timedelta(hours=16)
#                units =spill.units
#                skimmer1 = Skimmer(80, units=units, efficiency=0.36,
#                                  active_start=skim1_start,
#                                  active_stop=skim1_start + timedelta(hours=8))
#                skimmer2 = Skimmer(120, units=units, efficiency=0.2,
#                                  active_start=skim2_start,
#                                  active_stop=skim2_start + timedelta(hours=12))
#            #
#                chem_start = start_time + timedelta(hours=24)
#                c_disp = ChemicalDispersion(0.5, efficiency=0.4,
#                                            active_start=chem_start,
#                                            active_stop=chem_start + timedelta(hours=8),waves=waves)
#                burn_start = start_time + timedelta(hours=36)
#                burn = Burn(1000., .1, active_start=burn_start, efficiency=.2)
                model.environment += [water, wind,  waves]
#                model.weatherers += skimmer1
#                model.weatherers +=skimmer2
#                model.weatherers += burn
#                model.weatherers += c_disp
#            #    
                model.weatherers += Evaporation(water,wind)
                model.weatherers += Emulsification(waves)
                model.weatherers += NaturalDispersion(waves,water)
                
                return model
                
            if __name__ == '__main__':   
                model = make_model()
                for step in model:
                  print "step: %.4i -- memuse: %fMB" % (step['step_num'],
                                                   utilities.get_mem_use())
                                                   
            print '\nIteration number %i Completed'%q, 'for start time %i'%i, 'of Spill number%i'%s

#
print '\nA total of %i simulations performed'%number_of_simulations
print "It took %s to run the simulations" % (datetime.now() - simulation_start_time)

dir0='my_run_performance'
if not os.path.exists(dir0):
    os.makedirs(dir0)
my_run_performance = open(dir0+'/run_performance.txt','w')
my_run_performance.write("\nA total of %i simulations peformed" %number_of_simulations)
my_run_performance.write("\nIt took %s to run the simulations" % (datetime.now() - simulation_start_time))
my_run_performance.close()
#----------------------------------------------------------------------------------------
for i in range(1,3):
    print'\n','*******************'*i, 'PLEASE SCROLL'
    time.sleep(0.1)
#----------------------------------------------------------------------------------------
#    
    
print '\nNow Starting post-processing for CERTAIN TRAJECTORY'
#
import numpy
from netCDF4 import Dataset
from shapely.geometry import Point, Polygon

print "\nCalculation of MCF & Risk Matrix Is In Progress"

MCF=numpy.zeros((number_of_spills,number_of_starttimes,number_iterations_for_starttime, number_of_polygons))
risk_matrix=numpy.zeros((number_of_spills,number_of_starttimes,number_iterations_for_starttime, number_of_polygons))
processed_risk_matrix=[[[0.0 for k in xrange(8)] for y in xrange(number_of_polygons)] for x in xrange(number_of_spills)]
#

dir1='Text_results\CERTAIN/Details/Risk'
if not os.path.exists(dir1):
    os.makedirs(dir1)
    
dir2='Text_results\CERTAIN/Details/Processed_risk'
if not os.path.exists(dir2):
    os.makedirs(dir2)

dir3='Text_results\CERTAIN/Summary'
if not os.path.exists(dir3):
    os.makedirs(dir3)

dir4='Text_results\CERTAIN/Extreme'
if not os.path.exists(dir4):
    os.makedirs(dir4)

    
g = open(dir1+'/README.txt','w')
g.write("Row index is start-time number and column index is iteration number")
g.write("\n-----------------------------------------------------------------")
g.close()

q = open(dir2+'/README.txt','w')
q.write("From the first row to 8th row: MAX, MIN, PERCENTILE_10, PERCENTILE_50 ,PERCENTILE_90, AVERAGE, STANDARD DEVIATION, RELATIVE STANDARD DEVIATION")
q.write("\n--------------------------------------------------------------------------------------------------------------------------------------------")
q.close()

r = open(dir4+'/README.txt','w')
r.write("This file contains the risks that are overally greater than 90 precent of other risks; and the related indices:Overal 90th percentile, risk value,[spill_number,starttime_number,iteration_number,polygon_number]")
r.write("\n----------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
r.close()

s=open(dir4+'/Extreme_report.txt','w')


for spill_number in range (0,number_of_spills ):
    for starttime_number in range(0,number_of_starttimes):
        for iteration_number in range (0,number_iterations_for_starttime):
            post_nc=Dataset('Results\Spill %i\NC\Start_time %i\gulf_output%i.nc'%(spill_number, starttime_number,iteration_number),'r')
            lat=post_nc.variables['latitude']
            lon=post_nc.variables['longitude']
            mass=post_nc.variables['mass']
            id_number=post_nc.variables['id']
            particles_in_step=post_nc.variables['particle_count']
            # now initiating parameters
            steps=NC_steps+2
            for polygon_number in range(0,number_of_polygons):
                gulf_polygon=Polygon(P[polygon_number])
                j=2 # j=1 is for particles in the moment of release exactly in release location. not needed according to averaging time definition
                mass_particles=[] #mass
                c=[]
                while j<steps:
#                    print j
                    m=0
                    if j==2:lower_bound=0
                    else:lower_bound=numpy.sum(particles_in_step[0:j-2])
                    for i in range(lower_bound,(lower_bound+particles_in_step[j-1])):
                        
#                        print i,'(',lon[i] , lat[i], mass[i], ')'
                        point=Point (lon[i],lat[i])
                        if  gulf_polygon.contains(point):
                            m+=mass[i]     
                    mass_particles.append(m)
                    c.append(m/gulf_polygon.area)
                    j+=1
                MCF[spill_number,starttime_number,iteration_number,polygon_number]=numpy.average(c[:]) 
                risk_matrix[spill_number,starttime_number,iteration_number,polygon_number]=gulf_frequency_matrix[spill_number][starttime_number][iteration_number]*sensitivity[polygon_number]* MCF[spill_number,starttime_number,iteration_number,polygon_number]
#                print spill_number, starttime_number, iteration_number, polygon_number, risk_matrix[spill_number,starttime_number,iteration_number,polygon_number]
                print spill_number, starttime_number, iteration_number, polygon_number, datetime.now()

            
            post_nc.close()



##-----------------------------------------------------------------------------------------
#
#


    print "\n\nThe summary of information for polygons in Spill #%i:"%spill_number    
    for i in range (0,number_of_polygons):
        print "\nThe output summary for polygon#%i:"%i
        risk_max=numpy.max(risk_matrix[spill_number,:,:,i])
        print "\nMaximum risk:%f"%risk_max
        processed_risk_matrix[spill_number][i][0]=risk_max
        risk_min=numpy.min(risk_matrix[spill_number,:,:,i])
        print "Minimum risk:%f"%risk_min
        processed_risk_matrix[spill_number][i][1]=risk_min
        percentile_10=numpy.percentile(risk_matrix[spill_number,:,:,i],10)
        print "The 10th percentile:%f"%percentile_10
        processed_risk_matrix[spill_number][i][2]=percentile_10
        percentile_50=numpy.percentile(risk_matrix[spill_number,:,:,i],50)
        print "The 50th percentile (median):%f"%percentile_50
        processed_risk_matrix[spill_number][i][3]=percentile_50
        percentile_90=numpy.percentile(risk_matrix[spill_number,:,:,i],90)
        print "The 90th percentile:%f"%percentile_90
        processed_risk_matrix[spill_number][i][4]=percentile_90
        risk_mean=numpy.average(risk_matrix[spill_number,:,:,i])
        print "The average risk:%f"%risk_mean
        processed_risk_matrix[spill_number][i][5]=risk_mean
        risk_std=numpy.std(risk_matrix[spill_number,:,:,i])
        print "The standard deviation of risk:%f"%risk_std
        processed_risk_matrix[spill_number][i][6]=risk_std
        relative_std_percent=(risk_std/risk_mean)*100
        print "The relative standard deviation:%f\n\n"%relative_std_percent
        processed_risk_matrix[spill_number][i][7]=relative_std_percent
        
        numpy.savetxt(dir1+'/risk_S%iR%i.txt'%(spill_number,i),risk_matrix[spill_number,:,:,i])
        numpy.savetxt(dir2+'/processed_risk_S%iR%i.txt'%(spill_number,i),processed_risk_matrix[spill_number][i][:])

# -----------------------------------------------------------------------------------#        


# find overall maximums and the related indices to discover worst case inputs:-------#
extreme=[]
OVERAL_percentile_90=numpy.percentile(risk_matrix[:,:,:,:],90)
for spill_number in range (0,number_of_spills ):
    for starttime_number in range(0,number_of_starttimes):
        for iteration_number in range (0,number_iterations_for_starttime):
            for polygon_number in range(0,number_of_polygons):
                if risk_matrix[spill_number,starttime_number,iteration_number,polygon_number]>OVERAL_percentile_90:
                    print OVERAL_percentile_90, risk_matrix[spill_number,starttime_number,iteration_number,polygon_number],[spill_number,starttime_number,iteration_number,polygon_number]
                    extreme.append([OVERAL_percentile_90, risk_matrix[spill_number,starttime_number,iteration_number,polygon_number],[spill_number,starttime_number,iteration_number,polygon_number]])
for i in range(len(extreme)):
    s.write("%s\n"%extreme[i])
s.close()

##
#------------------------------------------------------------------------------------#
##AGGREGATION for RECEPTORS

print "\n\nThe Final Aggregation Result for Risk Receptor Polygons:"
max_risk_receptor=numpy.zeros(number_of_polygons,dtype=float)
min_risk_receptor=numpy.zeros(number_of_polygons,dtype=float)
percentile_10_receptor=numpy.zeros(number_of_polygons,dtype=float)
percentile_50_receptor=numpy.zeros(number_of_polygons,dtype=float)
percentile_90_receptor=numpy.zeros(number_of_polygons,dtype=float)
mean_risk_receptor=numpy.zeros(number_of_polygons,dtype=float)
#std_risk=numpy.zeros(number_of_polygons)
#std_relative_risk=numpy.zeros(number_of_polygons)

f = open(dir3+'/Results_summary.txt','w')
f.write("THE SUMMARY OF AGGREGATIVE RISK FOR POLYGONS")
f.write("\n--------------------------------------------")

for n in range(0,number_of_polygons):
    for s in range (0,number_of_spills):
        max_risk_receptor[n]+=numpy.sum(processed_risk_matrix[s][n][0])
        min_risk_receptor[n]+=numpy.sum(processed_risk_matrix[s][n][1])
        percentile_10_receptor[n]+=numpy.sum(processed_risk_matrix[s][n][2])
        percentile_50_receptor[n]+=numpy.sum(processed_risk_matrix[s][n][3])
        percentile_90_receptor[n]+=numpy.sum(processed_risk_matrix[s][n][4])
        mean_risk_receptor[n]+=numpy.sum(processed_risk_matrix[s][n][5])
#        
        
    print "\nAggregative max risk for polygon %i:"%n, max_risk_receptor[n]
    print "Aggregative min risk for polygon %i:"%n, min_risk_receptor[n]
    print "Aggregative 10th percentile for polygon %i:"%n, percentile_10_receptor[n]
    print "Aggregative 50th percentile for polygon %i:"%n, percentile_50_receptor[n]
    print "Aggregative 90th percentile for polygon %i:"%n, percentile_90_receptor[n]
    print "Aggregative mean risk for polygon %i:"%n, mean_risk_receptor[n]
    print 'TOTAL STD for Polygon %i:?'%n
    print 'TOTAL Relative STD for Polygon %i:?'%n
    

    f.write("\nPOLYGON %i:"%n)
    f.write("\nAggregative max risk for polygon %i:%f"%(n, max_risk_receptor[n]))
    f.write("\nAggregative min risk for polygon %i:%f"%(n, min_risk_receptor[n]))
    f.write("\nAggregative 10th percentile for polygon %i:%f"%(n, percentile_10_receptor[n]))
    f.write("\nAggregative 50th percentile for polygon %i:%f"%(n, percentile_50_receptor[n]))
    f.write("\nAggregative 90th percentile for polygon %i:%f"%(n, percentile_90_receptor[n]))
    f.write("\nAggregative mean risk for polygon %i:%f"%(n, mean_risk_receptor[n]))
    
    
##AGGREGATION for Sources (Spill locations)
print "\n\nThe Final Aggregation Result for Risk Sources(Spill Locations):"
max_risk_source=numpy.zeros(number_of_spills,dtype=float)
min_risk_source=numpy.zeros(number_of_spills,dtype=float)
percentile_10_source=numpy.zeros(number_of_spills,dtype=float)
percentile_50_source=numpy.zeros(number_of_spills,dtype=float)
percentile_90_source=numpy.zeros(number_of_spills,dtype=float)
mean_risk_source=numpy.zeros(number_of_spills,dtype=float)
#std_risk=numpy.zeros(number_of_spills)
#std_relative_risk=numpy.zeros(number_of_spills)

f.write("\n\nTHE SUMMARY OF AGGREGATIVE RISK FROM SOURCES")
f.write("\n--------------------------------------------")

for s in range (0,number_of_spills):
    for n in range(0,number_of_polygons):
        max_risk_source[s]+=numpy.sum(processed_risk_matrix[s][n][0])
        min_risk_source[s]+=numpy.sum(processed_risk_matrix[s][n][1])
        percentile_10_source[s]+=numpy.sum(processed_risk_matrix[s][n][2])
        percentile_50_source[s]+=numpy.sum(processed_risk_matrix[s][n][3])
        percentile_90_source[s]+=numpy.sum(processed_risk_matrix[s][n][4])
        mean_risk_source[s]+=numpy.sum(processed_risk_matrix[s][n][5])
#        
    print "\nAggregative max risk from source %i:"%s, max_risk_source[s]
    print "Aggregative min risk from source %i:"%s, min_risk_source[s]
    print "Aggregative 10th percentile from source %i:"%s, percentile_10_source[s]
    print "Aggregative 50th percentile from source %i:"%s, percentile_50_source[s]
    print "Aggregative 90th percentile from source %i:"%s, percentile_90_source[s]
    print "Aggregative mean risk from source %i:"%s, mean_risk_source[s]
    print 'TOTAL STD for source %i:?'%s
    print 'TOTAL Relative STD for source %i:?'%s
#    

    f.write("\nSOURCE %i:"%s)
    f.write("\nAggregative max risk from source %i:%f"%(s, max_risk_source[s]))
    f.write("\nAggregative min risk from source %i:%f"%(s, min_risk_source[s]))
    f.write("\nAggregative 10th percentile from source %i:%f"%(s, percentile_10_source[s]))
    f.write("\nAggregative 50th percentile from source %i:%f"%(s, percentile_50_source[s]))
    f.write("\nAggregative 90th percentile from source %i:%f"%(s, percentile_90_source[s]))
    f.write("\nAggregative mean risk from source %i:%f"%(s, mean_risk_source[s]))    
    
    
f.close()
    
    
    
#------------------------------------------------------------------------------------#
for i in range(1,3):
    print'\n', '*******************'*i, 'PLEASE SCROLL'
    time.sleep(0.1)
#------------------------------------------------------------------------------------#
    
print '\nNow Starting post-processing for UNCERTAIN TRAJECTORY'


print "\nCalculation of MCF & Risk Matrix Is In Progress"

MCF_uncertain=numpy.zeros((number_of_spills,number_of_starttimes,number_iterations_for_starttime, number_of_polygons))
risk_matrix_uncertain=numpy.zeros((number_of_spills,number_of_starttimes,number_iterations_for_starttime, number_of_polygons))
processed_risk_matrix_uncertain=[[[0.0 for l in xrange(8)] for y in xrange(number_of_polygons)] for x in xrange(number_of_spills)]
#

dir1_uncertain='Text_results\UNCERTAIN/Details/Risk'
if not os.path.exists(dir1_uncertain):
    os.makedirs(dir1_uncertain)
    
dir2_uncertain='Text_results\UNCERTAIN/Details/Processed_risk'
if not os.path.exists(dir2_uncertain):
    os.makedirs(dir2_uncertain)

dir3_uncertain='Text_results\UNCERTAIN/Summary'
if not os.path.exists(dir3_uncertain):
    os.makedirs(dir3_uncertain)
    
dir4_uncertain='Text_results\UNCERTAIN/Extreme'
if not os.path.exists(dir4_uncertain):
    os.makedirs(dir4_uncertain)
#
#    
g_uncertain = open(dir1_uncertain+'/README.txt','w')
g_uncertain.write("Row index is start-time number and column index is iteration number")
g_uncertain.write("\n-----------------------------------------------------------------")
g_uncertain.close()

q_uncertain = open(dir2_uncertain+'/README.txt','w')
q_uncertain.write("From the first row to 8th row: MAX, MIN, PERCENTILE_10, PERCENTILE_50 ,PERCENTILE_90, AVERAGE, STANDARD DEVIATION, RELATIVE STANDARD DEVIATION")
q_uncertain.write("\n--------------------------------------------------------------------------------------------------------------------------------------------")
q_uncertain.close()
    
r_uncertain= open(dir4_uncertain+'/README.txt','w')
r_uncertain.write("This file contains the risks that are overally greater than 90 precent of other risks; and the related indices:Overal 90th percentile, risk value,[spill_number,starttime_number,iteration_number,polygon_number]")
r_uncertain.write("\n----------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
r_uncertain.close()

s_uncertain=open(dir4_uncertain+'/Extreme_report.txt','w')


for spill_number in range (0,number_of_spills ):
    for starttime_number in range(0,number_of_starttimes):
        for iteration_number in range (0,number_iterations_for_starttime):
            post_nc_uncertain=Dataset('Results\Spill %i\NC\Start_time %i\gulf_output%i_uncertain.nc'%(spill_number, starttime_number,iteration_number),'r')
            lat=post_nc_uncertain.variables['latitude']
            lon=post_nc_uncertain.variables['longitude']
            mass=post_nc_uncertain.variables['mass']
            id_number=post_nc_uncertain.variables['id']
            particles_in_step=post_nc_uncertain.variables['particle_count']

            # now initiating parameters
            steps=NC_steps+2
            for polygon_number in range(0,number_of_polygons):
                gulf_polygon=Polygon(P[polygon_number])
                j=2
                mass_particles=[] #mass
                c=[]
                while j<steps:
                    m=0
                    if j==2:lower_bound=0
                    else:lower_bound=numpy.sum(particles_in_step[0:j-2])
                    for i in range(lower_bound,(lower_bound+particles_in_step[j-1])):
#                        print i,'(',lon[i] , lat[i], mass[i], ')'
                        point=Point (lon[i],lat[i])
                        if  gulf_polygon.contains(point):
                            m+=mass[i]     
                    mass_particles.append(m)
                    c.append(m/gulf_polygon.area)
                    j+=1
                MCF_uncertain[spill_number,starttime_number,iteration_number,polygon_number]=numpy.average(c[:]) 
                risk_matrix_uncertain[spill_number,starttime_number,iteration_number,polygon_number]=gulf_frequency_matrix[spill_number][starttime_number][iteration_number]*sensitivity[polygon_number]* MCF_uncertain[spill_number,starttime_number,iteration_number,polygon_number]
                print "uncertain", spill_number, starttime_number, iteration_number, polygon_number, datetime.now()

            
            post_nc_uncertain.close()
#-----------------------------------------------------------------------------------------


    print "\n\nThe summary of information for polygons in Spill #%i:"%spill_number    
    for i in range (0,number_of_polygons):
        print "\nThe output summary for polygon#%i:"%i
        risk_max=numpy.max(risk_matrix_uncertain[spill_number,:,:,i])
        print "\nMaximum risk:%f"%risk_max
        processed_risk_matrix_uncertain[spill_number][i][0]=risk_max
        risk_min=numpy.min(risk_matrix_uncertain[spill_number,:,:,i])
        print "Minimum risk:%f"%risk_min
        processed_risk_matrix_uncertain[spill_number][i][1]=risk_min
        percentile_10=numpy.percentile(risk_matrix_uncertain[spill_number,:,:,i],10)
        print "The 10th percentile:%f"%percentile_10
        processed_risk_matrix_uncertain[spill_number][i][2]=percentile_10
        percentile_50=numpy.percentile(risk_matrix_uncertain[spill_number,:,:,i],50)
        print "The 50th percentile (median):%f"%percentile_50
        processed_risk_matrix_uncertain[spill_number][i][3]=percentile_50
        percentile_90=numpy.percentile(risk_matrix_uncertain[spill_number,:,:,i],90)
        print "The 90th percentile:%f"%percentile_90
        processed_risk_matrix_uncertain[spill_number][i][4]=percentile_90
        risk_mean=numpy.average(risk_matrix_uncertain[spill_number,:,:,i])
        print "The average risk:%f"%risk_mean
        processed_risk_matrix_uncertain[spill_number][i][5]=risk_mean
        risk_std=numpy.std(risk_matrix_uncertain[spill_number,:,:,i])
        print "The standard deviation of risk:%f"%risk_std
        processed_risk_matrix_uncertain[spill_number][i][6]=risk_std
        relative_std_percent=(risk_std/risk_mean)*100
        print "The relative standard deviation:%f\n\n"%relative_std_percent
        processed_risk_matrix_uncertain[spill_number][i][7]=relative_std_percent
        
        numpy.savetxt(dir1_uncertain+'/risk_S%iR%i.txt'%(spill_number,i),risk_matrix_uncertain[spill_number,:,:,i])
        numpy.savetxt(dir2_uncertain+'/processed_risk_S%iR%i.txt'%(spill_number,i),processed_risk_matrix_uncertain[spill_number][i][:])

# -----------------------------------------------------------------------------------#        


# find overall maximums and the related indices to discover worst case inputs:-------#
extreme_uncertain=[]
OVERAL_percentile_90_uncertain=numpy.percentile(risk_matrix_uncertain[:,:,:,:],90)
for spill_number in range (0,number_of_spills ):
    for starttime_number in range(0,number_of_starttimes):
        for iteration_number in range (0,number_iterations_for_starttime):
            for polygon_number in range(0,number_of_polygons):
                if risk_matrix_uncertain[spill_number,starttime_number,iteration_number,polygon_number]>OVERAL_percentile_90_uncertain:
                    print OVERAL_percentile_90_uncertain, risk_matrix_uncertain[spill_number,starttime_number,iteration_number,polygon_number],[spill_number,starttime_number,iteration_number,polygon_number]
                    extreme_uncertain.append([OVERAL_percentile_90_uncertain, risk_matrix_uncertain[spill_number,starttime_number,iteration_number,polygon_number],[spill_number,starttime_number,iteration_number,polygon_number]])
for i in range(len(extreme_uncertain)):
    s_uncertain.write("%s\n"%extreme_uncertain[i])
s_uncertain.close()


#
#
#------------------------------------------------------------------------------------#
##AGGREGATION for RECEPTORS

print "\n\nThe Final Aggregation Result for Risk Receptor Polygons:"
max_risk_receptor_uncertain=numpy.zeros(number_of_polygons,dtype=float)
min_risk_receptor_uncertain=numpy.zeros(number_of_polygons,dtype=float)
percentile_10_receptor_uncertain=numpy.zeros(number_of_polygons,dtype=float)
percentile_50_receptor_uncertain=numpy.zeros(number_of_polygons,dtype=float)
percentile_90_receptor_uncertain=numpy.zeros(number_of_polygons,dtype=float)
mean_risk_receptor_uncertain=numpy.zeros(number_of_polygons,dtype=float)
#std_risk=numpy.zeros(number_of_polygons)
#std_relative_risk=numpy.zeros(number_of_polygons)

f_uncertain = open(dir3_uncertain+'/Results_summary.txt','w')
f_uncertain.write("THE SUMMARY OF AGGREGATIVE RISK FOR POLYGONS")
f_uncertain.write("\n--------------------------------------------")

for n in range(0,number_of_polygons):
    for s in range (0,number_of_spills):
        max_risk_receptor_uncertain[n]+=numpy.sum(processed_risk_matrix_uncertain[s][n][0])
        min_risk_receptor_uncertain[n]+=numpy.sum(processed_risk_matrix_uncertain[s][n][1])
        percentile_10_receptor_uncertain[n]+=numpy.sum(processed_risk_matrix_uncertain[s][n][2])
        percentile_50_receptor_uncertain[n]+=numpy.sum(processed_risk_matrix_uncertain[s][n][3])
        percentile_90_receptor_uncertain[n]+=numpy.sum(processed_risk_matrix_uncertain[s][n][4])
        mean_risk_receptor_uncertain[n]+=numpy.sum(processed_risk_matrix_uncertain[s][n][5])
#        
        
    print "\nAggregative max risk for polygon %i:"%n, max_risk_receptor_uncertain[n]
    print "Aggregative min risk for polygon %i:"%n, min_risk_receptor_uncertain[n]
    print "Aggregative 10th percentile for polygon %i:"%n, percentile_10_receptor_uncertain[n]
    print "Aggregative 50th percentile for polygon %i:"%n, percentile_50_receptor_uncertain[n]
    print "Aggregative 90th percentile for polygon %i:"%n, percentile_90_receptor_uncertain[n]
    print "Aggregative mean risk for polygon %i:"%n, mean_risk_receptor_uncertain[n]
    print 'TOTAL STD for Polygon %i:?'%n
    print 'TOTAL Relative STD for Polygon %i:?'%n
    
    f_uncertain.write("\nPOLYGON %i:"%n)
    f_uncertain.write("\nAggregative max risk for polygon %i:%f"%(n, max_risk_receptor_uncertain[n]))
    f_uncertain.write("\nAggregative min risk for polygon %i:%f"%(n, min_risk_receptor_uncertain[n]))
    f_uncertain.write("\nAggregative 10th percentile for polygon %i:%f"%(n, percentile_10_receptor_uncertain[n]))
    f_uncertain.write("\nAggregative 50th percentile for polygon %i:%f"%(n, percentile_50_receptor_uncertain[n]))
    f_uncertain.write("\nAggregative 90th percentile for polygon %i:%f"%(n, percentile_90_receptor_uncertain[n]))
    f_uncertain.write("\nAggregative mean risk for polygon %i:%f"%(n, mean_risk_receptor_uncertain[n]))    
    
f_uncertain.write("\n\nTHE SUMMARY OF AGGREGATIVE RISK FOR SOURCES")
f_uncertain.write("\n-------------------------------------------")
    
##AGGREGATION for Sources (Spill locations)
print "\n\nThe Final Aggregation Result for Risk Sources(Spill Locations):"
max_risk_source_uncertain=numpy.zeros(number_of_spills,dtype=float)
min_risk_source_uncertain=numpy.zeros(number_of_spills,dtype=float)
percentile_10_source_uncertain=numpy.zeros(number_of_spills,dtype=float)
percentile_50_source_uncertain=numpy.zeros(number_of_spills,dtype=float)
percentile_90_source_uncertain=numpy.zeros(number_of_spills,dtype=float)
mean_risk_source_uncertain=numpy.zeros(number_of_spills,dtype=float)
#std_risk=numpy.zeros(number_of_spills)
#std_relative_risk=numpy.zeros(number_of_spills)

for s in range (0,number_of_spills):
    for n in range(0,number_of_polygons):
        max_risk_source_uncertain[s]+=numpy.sum(processed_risk_matrix_uncertain[s][n][0])
        min_risk_source_uncertain[s]+=numpy.sum(processed_risk_matrix_uncertain[s][n][1])
        percentile_10_source_uncertain[s]+=numpy.sum(processed_risk_matrix_uncertain[s][n][2])
        percentile_50_source_uncertain[s]+=numpy.sum(processed_risk_matrix_uncertain[s][n][3])
        percentile_90_source_uncertain[s]+=numpy.sum(processed_risk_matrix_uncertain[s][n][4])
        mean_risk_source_uncertain[s]+=numpy.sum(processed_risk_matrix_uncertain[s][n][5])
#        
    print "\nAggregative max risk from source %i:"%s, max_risk_source_uncertain[s]
    print "Aggregative min risk from source %i:"%s, min_risk_source_uncertain[s]
    print "Aggregative 10th percentile from source %i:"%s, percentile_10_source_uncertain[s]
    print "Aggregative 50th percentile from source %i:"%s, percentile_50_source_uncertain[s]
    print "Aggregative 90th percentile from source %i:"%s, percentile_90_source_uncertain[s]
    print "Aggregative mean risk from source %i:"%s, mean_risk_source_uncertain[s]
    print 'TOTAL STD for source %i:?'%s
    print 'TOTAL Relative STD for source %i:?'%s

    f_uncertain.write("\nSOURCE %i:"%s)
    f_uncertain.write("\nAggregative max risk from source %i:%f"%(s, max_risk_source_uncertain[s]))
    f_uncertain.write("\nAggregative min risk from source %i:%f"%(s, min_risk_source_uncertain[s]))
    f_uncertain.write("\nAggregative 10th percentile from source %i:%f"%(s, percentile_10_source_uncertain[s]))
    f_uncertain.write("\nAggregative 50th percentile from source %i:%f"%(s, percentile_50_source_uncertain[s]))
    f_uncertain.write("\nAggregative 90th percentile from source %i:%f"%(s, percentile_90_source_uncertain[s]))
    f_uncertain.write("\nAggregative mean risk from source %i:%f"%(s, mean_risk_source_uncertain[s]))
    
f_uncertain.close()


#----------------------------------------------------------------------------------------------------------------------

print "\nIt took %s to run the simulations and to save & process the results" % (datetime.now() - TOTAL_start_time)

my_run_performance = open(dir0+'/run_performance.txt','a+')
my_run_performance.write("\n\nIt took %s to run the simulations and to save & process the results" % (datetime.now() - TOTAL_start_time))
my_run_performance.write("\n-----------------------------------------------------------------")
my_run_performance.close()


# -----------------------------------------------------------------------------------#
#Visualization
print "\n\nRisk histogram:"
plt.hist(risk_matrix[0,:,:,0])

#------------------------------------------------------------------------------------#

