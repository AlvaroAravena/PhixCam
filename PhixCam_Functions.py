import numpy as np
import matplotlib.pyplot as plt
import os
import elevation
import tifffile
import utm
from scipy import interpolate
from matplotlib import transforms
import cv2
from tkinter import filedialog
from math import sin , cos , sqrt , atan2 , radians , log , factorial , tan
import sys
from tkinter import messagebox

def interpol_pos_max( lon1 , lat1 , step_lon_deg , step_lat_deg , lon_cen , lat_cen , cells_lon , cells_lat , Topography ):

	dlon = int( np.floor( ( lon_cen - lon1 ) / step_lon_deg ) )
	dlat = ( cells_lat - 2 ) - int( np.floor( ( lat_cen - lat1 ) / step_lat_deg ) )
	if( dlon >= ( cells_lon - 1.0 ) or dlat >= ( cells_lat - 1.0 ) or dlon < 0.0 or dlat < 0.0 ):
		return -9999
	aux_lon = 2.0 * ( lon_cen - ( dlon * step_lon_deg + lon1 ) - step_lon_deg / 2.0 ) / step_lon_deg
	aux_lat = 2.0 * ( - lat_cen + ( ( cells_lat - 1.0 - dlat ) * step_lat_deg + lat1 ) - step_lat_deg / 2.0 ) / step_lat_deg
	dc = np.max( [ Topography[ dlat ][ dlon ] , Topography[ dlat ][ dlon + 1 ] , Topography[ dlat + 1 ][ dlon ] , Topography[ dlat + 1 ][ dlon + 1 ] ] )
	return ( dc )

def interpol_pos( lon1 , lat1 , step_lon_deg , step_lat_deg , lon_cen , lat_cen , cells_lon , cells_lat , Topography ):

	dlon = int( np.floor( ( lon_cen - lon1 ) / step_lon_deg ) )
	dlat = ( cells_lat - 2 ) - int( np.floor( ( lat_cen - lat1 ) / step_lat_deg ) )
	if( dlon >= ( cells_lon - 1.0 ) or dlat >= ( cells_lat - 1.0 ) or dlon < 0.0 or dlat < 0.0 ):
		return -9999
	aux_lon = 2.0 * ( lon_cen - ( dlon * step_lon_deg + lon1 ) - step_lon_deg / 2.0 ) / step_lon_deg
	aux_lat = 2.0 * ( - lat_cen + ( ( cells_lat - 1.0 - dlat ) * step_lat_deg + lat1 ) - step_lat_deg / 2.0 ) / step_lat_deg
	dc = ( Topography[ dlat ][ dlon ] + Topography[ dlat ][ dlon + 1 ] + Topography[ dlat + 1 ][ dlon ] + Topography[ dlat + 1 ][ dlon + 1 ] ) / 4
	[ x3 , y3 , z3 ] = [ 0.0 , 0.0 , dc ]
	if( aux_lon >= 0.0 and abs( aux_lon ) >= abs( aux_lat ) ):
		[ x1 , y1 , z1 ] = [ 1.0 , 1.0 , Topography[ dlat + 1 ][ dlon + 1 ] ] 
		[ x2 , y2 , z2 ] = [ 1.0 , -1.0 , Topography[ dlat ][ dlon + 1 ] ] 
	elif( aux_lat >= 0.0 and abs( aux_lon ) < abs( aux_lat ) ):
		[ x1 , y1 , z1] = [ -1.0 , 1.0 , Topography[ dlat + 1 ][ dlon ] ] 
		[ x2 , y2 , z2] = [ 1.0 , 1.0 , Topography[ dlat + 1 ][ dlon + 1 ] ] 
	elif( aux_lon < 0.0 and abs( aux_lon ) >= abs( aux_lat ) ):
		[ x1 , y1 , z1 ] = [ -1.0 , 1.0 , Topography[ dlat + 1 ][ dlon ] ] 
		[ x2 , y2 , z2 ] = [ -1.0 , -1.0 , Topography[ dlat ][ dlon ] ] 
	else:
		[ x1 , y1 , z1 ] = [ -1.0 , -1.0 , Topography[ dlat ][ dlon ] ] 
		[ x2 , y2 , z2 ] = [ 1.0 , -1.0 , Topography[ dlat ][ dlon + 1 ] ]
	f1 = ( y2 - y1 ) * ( z3 - z1 ) - ( y3 - y1 ) * ( z2 - z1 )
	f2 = ( z2 - z1 ) * ( x3 - x1 ) - ( z3 - z1 ) * ( x2 - x1 )
	f3 = ( x2 - x1 ) * ( y3 - y1 ) - ( x3 - x1 ) * ( y2 - y1 )

	return ( ( - aux_lon * f1 - aux_lat * f2 ) / f3 + dc )

def import_map( toponame , lon , lat ):
	cam = [ lon , lat ]
	exe_path = os.getcwd()
	dist_range = 160000.0
	r_lat = 0.5 * dist_range / distance_two_points( cam[1] - 0.5 , cam[1] + 0.5 , cam[0] , cam[0] )
	r_lon = 0.5 * dist_range / distance_two_points( cam[1] , cam[1] , cam[0] - 0.5 , cam[0] + 0.5 )
	lon1 = np.round( cam[0] - r_lon , 1 )
	lon2 = np.round( cam[0] + r_lon , 1 )
	lat1 = np.round( cam[1] - r_lat , 1 )
	lat2 = np.round( cam[1] + r_lat , 1 )
	aux_lon = np.array( [ lon1 , lon2 ] )
	aux_lat = np.array( [ lat1 , lat2 ] )
	lon1 = min( aux_lon )
	lon2 = max( aux_lon )
	lat1 = min( aux_lat )
	lat2 = max( aux_lat )
	file_txt = open( os.path.join( exe_path , 'Cities.txt' ) )
	line = file_txt.readlines()
	file_txt.close()
	for population in range( 10000 , 10000000 , 10000 ):
		Cities = []
		for i in range( 1 , len( line ) ):
			aux = line[ i ].split( ',' )
			pop = float( aux[ 4 ] )
			lat_dat = float( aux[ 5 ] )
			lon_dat = float( aux[ 6 ] )
			if( lon_dat > lon1 and lon_dat < lon2 and lat_dat > lat1 and lat_dat < lat2 and pop > population ):
				Cities.append( [ lon_dat , lat_dat , aux[ 2 ] ] )
		if( len( Cities ) <= 5 ):
			break
	elevation.clip( bounds = ( lon1 , lat1 , lon2 , lat2 ) , output = exe_path + '/' + toponame + '.tif' )
	fp = toponame + '.tif'
	image = tifffile.imread( fp )
	elevation.clean()
	Topography = np.array( image )
	Topography_Sea = Topography + 0.0
	Topography_Sea[ Topography_Sea[ : , : ] <= 0 ] = -1.0 * np.sqrt( -1.0 * Topography_Sea[ Topography_Sea[ : , : ] <= 0 ] )
	Topography_Sea[ Topography_Sea[ : , : ] > 0 ] = np.nan
	Topography_Sea = Topography_Sea * -1.0
	Topography = ( Topography + abs( Topography ) ) / 2.0
	cells_lon = Topography.shape[ 1 ]
	cells_lat = Topography.shape[ 0 ]
	text_file = open( 'Topographies/' + toponame + '.txt' , 'w' )
	text_file.write( 'lon1 ' + str( lon1 ) + '\n' )
	text_file.write( 'lon2 ' + str( lon2 ) + '\n' )
	text_file.write( 'lat1 ' + str( lat1 ) + '\n' )
	text_file.write( 'lat2 ' + str( lat2 ) + '\n' )
	text_file.write( 'cells_lon ' + str( cells_lon ) + '\n' )
	text_file.write( 'cells_lat ' + str( cells_lat ) + '\n' )
	for i in range( cells_lat ):
		for j in range( cells_lon ):
			text_file.write( str( Topography[ i , j ] ) + ' ' )
		text_file.write( '\n' )
	text_file.close()
	return [ lon1 , lon2 , lat1 , lat2 , Cities , Topography , Topography_Sea , cells_lon , cells_lat ]

def read_map_deg( current_path , topography_file ):
	exe_path = os.getcwd()
	try:
 		file_txt = open( topography_file )
	except:
 		print( topography_file + ' not found in ' + str( current_path ) + '.' )
 		sys.exit( 0 )
	line = file_txt.readlines()
	file_txt.close()
	[ lon1 , lon2 , lat1 , lat2 , cells_lon , cells_lat ] = [ np.nan , np.nan , np.nan , np.nan , np.nan , np.nan ]
	for i in range( 0 , 10 ):
		aux = line[ i ].split()
		if( aux[ 0 ] == 'lon1' ):
			lon1 = float( aux[ 1 ] )
		if( aux[ 0 ] == 'lon2' ):
			lon2 = float( aux[ 1 ] )
		if( aux[ 0 ] == 'lat1' ):
			lat1 = float( aux[ 1 ] )
		if( aux[ 0 ] == 'lat2' ):
			lat2 = float( aux[ 1 ] )
		if( aux[ 0 ] == 'cells_lon' ):
			cells_lon = int( aux[ 1 ] )
		if( aux[ 0 ] == 'cells_lat' ):
			cells_lat = int( aux[ 1 ] )
		if( len( aux ) >= 10 ):
			indexini = i
			break
	if( np.isnan( lon1 ) or np.isnan( lon2 ) or np.isnan( lat1 ) or np.isnan( lat2 ) or np.isnan( cells_lon ) or np.isnan( cells_lat ) ):
		print( 'Problems in topography file.' )
		sys.exit( 0 )
	Topography = np.zeros( ( cells_lat , cells_lon ) )
	for i in range( indexini , indexini + cells_lat ):
		aux = line[ i ].split()
		for j in range( 0 , cells_lon ):
			Topography[ i - indexini , j ] = float( aux[ j ] )
	Topography_Sea = Topography + 0.0
	Topography_Sea[ Topography_Sea[ : , : ] <= 0] = -1.0 * np.sqrt( -1.0 * Topography_Sea[ Topography_Sea[ : , : ] <= 0 ] )
	Topography_Sea[ Topography_Sea[ : , : ] > 0] = np.nan
	Topography_Sea = Topography_Sea * -1.0
	Topography = ( Topography + abs( Topography ) ) / 2.0
	file_txt = open( os.path.join( exe_path , 'Cities.txt' ) )
	line = file_txt.readlines()
	file_txt.close()
	for population in range( 10000 , 10000000 , 10000 ):
		Cities = []
		for i in range( 1 , len( line ) ):
			aux = line[ i ].split( ',' )
			pop = float( aux[ 4 ] )
			lat_dat = float( aux[ 5 ] )
			lon_dat = float( aux[ 6 ] )
			if( lon_dat > lon1 and lon_dat < lon2 and lat_dat > lat1 and lat_dat < lat2 and pop > population ):
				Cities.append( [ lon_dat , lat_dat , aux[ 2 ] ] )
		if( len( Cities ) <= 5 ):
			break
	return [ lon1 , lon2 , lat1 , lat2 , Cities , Topography , Topography_Sea , cells_lon , cells_lat ]

def read_map_utm( current_path , topography_file ):
	try:
		file_txt = open( topography_file )
	except:
		print( topography_file + ' not found in ' + str( current_path ) + '.' )
		sys.exit( 0 )
	line = file_txt.readlines()
	file_txt.close()
	[ n_north , n_east , cellsize , east_cor , north_cor , nodata_value ] = [ np.nan , np.nan , np.nan , np.nan , np.nan , np.nan ]
	for i in range( 0 , 10 ):
		aux = line[ i ].split()
		if( aux[ 0 ] == 'nrows' ):
			n_north = int( aux[ 1 ] )
		if( aux[ 0 ] == 'ncols' ):
			n_east = int( aux[ 1 ] )
		if( aux[ 0 ] == 'cellsize' ):
			cellsize = float( aux[ 1 ] )
		if( aux[ 0 ] == 'xllcorner' ):
			east_cor = float( aux[ 1 ] )
		if( aux[ 0 ] == 'yllcorner' ):
			north_cor = float( aux[ 1 ] )
		if( aux[ 0 ] == 'NODATA_value' ):
			nodata_value = float( aux[ 1 ] )
		if( len( aux ) >= 10 ):
			indexini = i
			break
	if( np.isnan( n_north ) or np.isnan( n_east ) or np.isnan( cellsize ) or np.isnan( east_cor ) or np.isnan( north_cor ) ):
		print( 'Problems in topography file.' )
		sys.exit( 0 )
	Topography = np.zeros( ( n_north , n_east ) )
	for i in range( indexini , indexini + n_north ):
		aux = line[ i ].split()
		for j in range( 0 , n_east ):
			Topography[ i - indexini , j ] = float( aux[ j ] )
	if( not np.isnan( nodata_value ) ):
		aux = np.where( Topography == nodata_value )
		Topography[ np.where( Topography == nodata_value ) ] = np.nan
	Topography_Sea = Topography + 0.0
	Topography_Sea[ Topography_Sea[ : , : ] <= 0 ] = -1.0 * np.sqrt( -1.0 * Topography_Sea[ Topography_Sea[ : , : ] <= 0 ] )
	Topography_Sea[ Topography_Sea[ : , : ] > 0 ] = np.nan
	Topography_Sea = Topography_Sea * -1.0
	Topography = ( Topography + abs( Topography ) ) / 2.0

	return [ Topography , Topography_Sea , n_north , n_east , cellsize , east_cor , north_cor ]
	
def plot_deg( lon , lat , lon1 , lon2 , lat1 , lat2 , Cities , Topography , Topography_Sea , cells_lon , cells_lat , ang1 , ang2 ):

	utm1 = utm.from_latlon( lat1 , lon1 )
	utm2 = utm.from_latlon( lat2 , lon2 )
	if( utm1[ 2 ] == utm2[ 2 ] and utm1[ 3 ] == utm2[ 3 ] ):
		distance_lon = abs( utm2[ 0 ] - utm1[ 0 ] )
		distance_lat = abs( utm2[ 1 ] - utm1[ 1 ] )
	else:
		distance_lon = distance_two_points( lat1 , lat1 , lon1 , lon2 )
		distance_lat = distance_two_points( lat1 , lat2 , lon1 , lon1 )
	utm_save = utm.from_latlon( min( lat1 , lat2 ) , min( lon1 , lon2 ) )
	step_lon_m = distance_lon / ( cells_lon - 1 )
	step_lat_m = distance_lat / ( cells_lat - 1 )
	matrix_lon = np.zeros( ( cells_lat , cells_lon ) )
	matrix_lat = np.zeros( ( cells_lat , cells_lon ) )
	for i in range( 0 , cells_lon ): 
		matrix_lon[ : , i ] = lon1 + ( lon2 - lon1 ) * ( i ) / ( cells_lon - 1 )
	for j in range( 0 , cells_lat ):
		matrix_lat[ j , : ] = lat1 + ( lat2 - lat1 ) * ( cells_lat - 1 - j ) / ( cells_lat - 1 )
	step_lon_deg = ( lon2 - lon1 ) / ( cells_lon - 1 )
	step_lat_deg = ( lat2 - lat1 ) / ( cells_lat - 1 )
	plt.figure( 1 , figsize = ( 8.0 , 5.0 ) )
	plt.axes().set_aspect( step_lat_m / step_lon_m )
	cmapg = plt.cm.get_cmap( 'Greys' )
	cmaps = plt.cm.get_cmap( 'Blues' ) 
	CS_Topo = plt.contourf( matrix_lon , matrix_lat , Topography , 100 , alpha = 1.0 , cmap = cmapg , antialiased = True )
	CS_Sea = plt.contourf( matrix_lon , matrix_lat , Topography_Sea , 100 , alpha = 0.5 , cmap = cmaps , antialiased = True )
	plt.xlabel( 'Longitude $[º]$' )
	plt.ylabel( 'Latitude $[º]$' )
	plt.xlim( lon1 , lon2 )
	plt.ylim( lat1 , lat2 )
	for i in range( len( Cities ) ):
		plt.text( float( Cities[ i ][ 0 ] ) , float( Cities[ i ][ 1 ] ) , str( Cities[ i ][ 2 ] ) , horizontalalignment = 'center' , verticalalignment ='center' , fontsize = 6 )
	plt.plot( lon , lat , 'b.' , markersize = 5 )
	if( ang2 - ang1 < 360.0 ):
		plt.plot( [ lon , lon + 0.5 * np.sin( ang1 * np.pi / 180) ] , [ lat , lat + 0.5 * np.cos( ang1 * np.pi / 180) ] , 'b:' )
		plt.plot( [ lon , lon + 0.5 * np.sin( ang2 * np.pi / 180) ] , [ lat , lat + 0.5 * np.cos( ang2 * np.pi / 180) ] , 'b:' )
	plt.show()

def plot_utm( east , north , east_cor , north_cor , n_east , n_north , cellsize , Topography , Topography_Sea , ang1 , ang2 ):

	matrix_north = np.zeros( ( n_north , n_east ) )
	matrix_east = np.zeros( ( n_north , n_east ) )
	for i in range( 0 , n_east ):
		matrix_east[ : , i ] = ( east_cor + cellsize * i )
	for j in range( 0 , n_north ):
		matrix_north[ j , : ] = ( north_cor + cellsize * j )
	matrix_north = matrix_north[ range( len( matrix_north[ : , 0 ] ) -1 , -1 , -1 ) , : ]
	plt.figure( 1 , figsize = ( 8.0 , 5.0 ) )
	plt.axes().set_aspect( 1.0 )
	cmapg = plt.cm.get_cmap( 'Greys' )
	cmaps = plt.cm.get_cmap( 'Blues' ) 
	CS_Topo = plt.contourf( matrix_east , matrix_north , Topography , 100 , alpha = 1.0 , cmap = cmapg , antialiased = True )
	CS_Sea = plt.contourf( matrix_east , matrix_north , Topography_Sea , 100 , alpha = 0.5 , cmap = cmaps , antialiased = True )
	plt.xlabel( 'East [m]' )
	plt.ylabel( 'North [m]' )
	plt.xlim( east_cor , east_cor + cellsize * ( n_east - 1 ) )
	plt.ylim( north_cor , north_cor + cellsize * ( n_north - 1 ) )
	plt.plot( east , north , 'b.' , markersize = 5 )
	if( ang2 - ang1 < 360.0 ):
		plt.plot( [ east , east + 10000 * np.sin( ang1 * np.pi / 180) ] , [ north , north + 10000 * np.cos( ang1 * np.pi / 180) ] , 'b:' )
		plt.plot( [ east , east + 10000 * np.sin( ang2 * np.pi / 180) ] , [ north , north + 10000 * np.cos( ang2 * np.pi / 180) ] , 'b:' )
	plt.show()
		
def normalize( datax , datay ):
	if( ( np.max( datax ) - np.min( datax ) ) > 0 ):
		datay = ( datay - np.min( datay ) ) / ( np.max( datax ) - np.min( datax ) )
		datax = ( datax - np.min( datax ) ) / ( np.max( datax ) - np.min( datax ) )
	return [ datax , datay ]

def rotate( inclination ):
	radinclination = np.pi * inclination / 180.0
	return np.array([[np.cos(radinclination), -np.sin(radinclination)], [np.sin(radinclination), np.cos(radinclination)]])
    
def normalize_inclined( datax , datay , ninclination , inclination_1 , inclination_2 ):
	inclination_steps = np.arange( inclination_1 , inclination_2 + ( inclination_2 - inclination_1 ) / ( ninclination - 1 ) / 2 , ( inclination_2 - inclination_1  ) / ( ninclination - 1 ) )
	[ x0 , y0 ] = [ datax[0] , datay[0] ]
	list_datax = []
	list_datay = []
	for i in inclination_steps:
		new_xy = np.column_stack( ( datax - x0 , datay - y0 )) @ rotate( i ) + [ x0 , y0 ]
		[ cur_datax , cur_datay ] = normalize( new_xy[:,0] , new_xy[:,1] )
		list_datax.append( cur_datax )
		list_datay.append( cur_datay )				
	return [ list_datax , list_datay ]

def normalize_inclined_list( datax , datay , list_inclination ):
	[ x0 , y0 ] = [ datax[0] , datay[0] ]
	list_datax = []
	list_datay = []
	for i in list_inclination:
		new_xy = np.column_stack( ( datax - x0 , datay - y0 )) @ rotate( i ) + [ x0 , y0 ]
		[ cur_datax , cur_datay ] = normalize( new_xy[:,0] , new_xy[:,1] )
		list_datax.append( cur_datax )
		list_datay.append( cur_datay )				
	return [ list_datax , list_datay ]

def candidates_initial( candidates_x , candidates_y , candidates_z , candidates_sens , inclination , stats_data , angleshor , anglesver , view_init , windows , sens , stepang ):
	limit_points = 0.1 * len( angleshor )
	for i in range( len( view_init ) ):
		for j in range( len( windows ) ):
			if( view_init[i] + windows[j] <= np.max( angleshor ) ):
				cur_indexes = np.where( np.abs( angleshor - view_init[i] - windows[j] / 2 ) <=  windows[j] / 2 )
				[x,y] = normalize( angleshor[ cur_indexes ] , anglesver[ cur_indexes ] )
				if( np.abs( np.max(y) - stats_data[0] ) <= sens and len( cur_indexes[ 0 ] ) > limit_points ):
					cur_d1 = np.polyfit( x, y, 1 )[0]
					if( np.abs( cur_d1 - stats_data[1] ) <= sens ):
						cur_d2 = np.polyfit( x, y, 2 )[0]
						if( np.abs( cur_d2 - stats_data[2] ) <= sens ):
							candidates_x = np.append( candidates_x , view_init[i] )
							candidates_y = np.append( candidates_y , windows[j] )
							candidates_z = np.append( candidates_z , inclination )
							candidates_sens = np.append( candidates_sens , 2 * max([ np.abs( np.max(y) - stats_data[0] ) , np.abs( cur_d1 - stats_data[1] ) , np.abs( cur_d2 - stats_data[2] ) ] ) )
							sens = min([ sens ,  max([ 0.1 , 2 * np.abs( np.max(y) - stats_data[0] ) , 2 * np.abs( cur_d1 - stats_data[1] ) , 2 * np.abs( cur_d2 - stats_data[2] ) ] ) ] )
	return [ candidates_x , candidates_y , candidates_z , candidates_sens , sens ]

def candidates_zoom( candidates_x , candidates_y , candidates_z , stepang , stepinclination , zoomfactor , angleshor , anglesver , stats_data , sens , bordex , bordey , error_zoom ):
	candidates_x_zoom = np.zeros( 0 )
	candidates_y_zoom = np.zeros( 0 )
	candidates_inclination_zoom  = np.zeros( 0 )
	limit_points = 50
	for i in range( len( candidates_x ) ):
		list_view_init = np.arange( candidates_x[ i ] - stepang , candidates_x[ i ] + stepang + stepang / zoomfactor / 2 , stepang / zoomfactor )
		list_windows = np.arange( candidates_y[ i ] - stepang , candidates_y[ i ] + stepang + stepang / zoomfactor / 2 , stepang / zoomfactor )
		list_inclinations = np.arange( candidates_z[ i ] - stepinclination , candidates_z[ i ] + stepinclination + stepinclination / zoomfactor / 2 , stepinclination / zoomfactor )
		[ list_datax , list_datay ] = normalize_inclined_list( bordex , -bordey , list_inclinations )
		for j in list_view_init:
			for k in list_windows:
				if( j + k <= np.max( angleshor ) ):
					cur_indexes = np.where( np.abs( angleshor - j - k / 2 ) <=  k / 2 )
					if( len( cur_indexes[0] ) > limit_points ):
						[x,y] = normalize( angleshor[ cur_indexes ] , anglesver[ cur_indexes ] )
						if( np.abs( np.max(y) - stats_data[0] ) <= sens ):
							cur_d1 = np.polyfit( x , y , 1 )[0]
							if( np.abs( cur_d1 - stats_data[1] ) <= sens ):
								cur_d2 = np.polyfit( x, y, 2 )[0]
								if( np.abs( cur_d2 - stats_data[2] ) <= sens ):
									f = interpolate.interp1d( x , y )
									for m in range( len( list_inclinations ) ):
										datay_new = f( list_datax[ m ] )
										shift = np.mean( list_datay[ m ] ) - np.mean( datay_new )
										datay_new = datay_new + shift
										sumerr = np.sum( ( list_datay[ m ] - datay_new ) * ( list_datay[ m ] - datay_new ) )
										if( error_zoom > sumerr ):
											error_zoom = sumerr
											candidates_x_zoom = np.append( candidates_x_zoom , j )
											candidates_y_zoom = np.append( candidates_y_zoom , k )
											candidates_inclination_zoom = np.append( candidates_inclination_zoom , list_inclinations[ m ] )
		if( i % 20 == 0 ):
			print('Second stage: ' + str( round( i / len( candidates_x ) * 100 , 2 ) ) + '% finished.' )
	return [ candidates_x_zoom , candidates_y_zoom , candidates_inclination_zoom , error_zoom ]

def get_profile_deg( run_name , points , lon , lat , lon1 , lon2 , lat1 , lat2 , Cities, Topography , Topography_Sea , cells_lon , cells_lat , ang1 , ang2 , h0 , step_distance , max_distance ):
	step_angles = np.abs( ang2 - ang1 ) / ( points )
	cam = [ lon , lat ]
	step_lon_deg = ( lon2 - lon1 ) / ( cells_lon - 1 )
	step_lat_deg = ( lat2 - lat1 ) / ( cells_lat - 1 )
	camutm = utm.from_latlon( cam[ 1 ] , cam[ 0 ] )
	angles = np.arange( - 1 * ang1 + 90.0 , - 1 * ang2 + 90.0 , -1.0 * step_angles )
	angles = angles * np.pi / 180.0
	verangles = np.arange( 0.0 , 90.0 , step_angles )
	verangles = verangles * np.pi / 180.0
	tanverangles = np.tan( verangles )
	verangles_2 = np.arange( 0.0 , 90.0 , -0.1 )
	verangles_2 = verangles_2 * np.pi / 180.0
	datview = np.zeros( len( angles ) )
	for i in range( len( angles ) ):
		d = np.arange( 0 , 1000 * max_distance , step_distance )
		h = np.arange( len(d) )
		current_utm = [ camutm[0] , camutm[1] ]
		current_deg = utm.to_latlon( current_utm[0] , current_utm[1] , camutm[2] , camutm[3] ) 
		h[ 0 ] = interpol_pos( lon1 , lat1 , step_lon_deg , step_lat_deg , current_deg[ 1 ] , current_deg[ 0 ] , cells_lon , cells_lat , Topography )
		if( h[ 0 ] == -9999 ):
		    h[ 0 ] = 0
		counter = 1 
		for distance in d[1:]:
			current_utm = [ camutm[0] + distance * np.cos( angles[i] ) , camutm[1] + distance * np.sin( angles[i] ) ]
			current_deg = utm.to_latlon( current_utm[0] , current_utm[1] , camutm[2] , camutm[3] ) 
			h[ counter ] = interpol_pos( lon1 , lat1 , step_lon_deg , step_lat_deg , current_deg[ 1 ] , current_deg[ 0 ] , cells_lon , cells_lat , Topography )
			counter += 1
		index_init = [ 0 ]
		for j in range( len( verangles_2 ) ):
			hv = h[0] + h0 + d * tanverangles_2[j]
			indexes = np.where( hv < h )
			if( indexes[0].size > 0 ):
				index_init = np.where( verangles >= verangles_2[j] )[0]
				break
		for j in range( index_init[0] , len( verangles ) ):
			hv = h[0] + h0 + d * tanverangles[j]
			indexes = np.where( hv < h )
			if( indexes[0].size > 0 ):
				datview[i] = verangles[j]
			else:
				break
		if( i % 100  == 0 ):
			print( "  Completed: " + str( np.round( 100 *  i / len( angles ) , 2 ) ) + "%." )
	angles_profile = np.append( lon , np.append( np.append( h[0] + h0 , h0 ) , np.arange( ang1 , ang2 , step_angles ) ) )
	datview = np.append( lat , np.append( np.append( 1 , ang2 ) , datview * 180 / np.pi ) )
	if( run_name != '' ):
		np.savetxt( 'Horizons/' + run_name + '.txt' , np.transpose( [ angles_profile , datview ] ) , fmt='%.6f' ) 
	return [ angles_profile[3:] , datview [3:] , h[0] + h0 ]

def get_profile_utm( run_name , points , east , north , east_cor , north_cor , cellsize , Topography , n_east , n_north , ang1 , ang2 , h0 , step_distance , max_distance ):
	step_angles = np.abs( ang2 - ang1 ) / ( points )
	camutm = [ east , north ]
	angles = np.arange( - 1 * ang1 + 90.0 , - 1 * ang2 + 90.0 , -1.0 * step_angles )
	angles = angles * np.pi / 180.0
	verangles = np.arange( 0.0 , 90.0 , step_angles )
	verangles = verangles * np.pi / 180.0
	tanverangles = np.tan( verangles )
	verangles_2 = np.arange( 0.0 , 90.0 , -0.1 )
	verangles_2 = verangles_2 * np.pi / 180.0
	datview = np.zeros( len( angles ) )
	for i in range( len( angles ) ):
		d = np.arange( 0 , 1000 * max_distance , step_distance )
		h = np.arange( len(d) )
		current_utm = [ camutm[0] , camutm[1] ]
		h[ 0 ] = interpol_pos( east_cor , north_cor , cellsize , cellsize , current_utm[ 0 ] , current_utm[ 1 ] , n_east , n_north , Topography )
		if( h[ 0 ] == -9999 ):
		    h[ 0 ] = 0
		counter = 1 
		for distance in d[1:]:
			current_utm = [ camutm[0] + distance * np.cos( angles[i] ) , camutm[1] + distance * np.sin( angles[i] ) ]
			h[ counter ] = interpol_pos( east_cor , north_cor , cellsize , cellsize , current_utm[ 0 ] , current_utm[ 1 ] , n_east , n_north , Topography )
			counter += 1
		index_init = [ 0 ]
		for j in range( len( verangles_2 ) ):
			hv = h[0] + h0 + d * tanverangles_2[j]
			indexes = np.where( hv < h )
			if( indexes[0].size > 0 ):
				index_init = np.where( verangles >= verangles_2[j] )[0]
				break
		for j in range( index_init[0] , len( verangles ) ):
			hv = h[0] + h0 + d * tanverangles[j]
			indexes = np.where( hv < h )
			if( indexes[0].size > 0 ):
				datview[i] = verangles[j]
			else:
				break
		if( i % 100  == 0 ):
			print( "  Completed: " + str( np.round( 100 *  i / len( angles ) , 2 ) ) + "%." )
	angles_profile = np.append( east , np.append( np.append( h[0] + h0 , h0 ) , np.arange( ang1 , ang2 , step_angles ) ) )
	datview = np.append( north , np.append( np.append( 2 , ang2 ) , datview * 180 / np.pi ) )
	if( run_name != '' ):
		np.savetxt( 'Horizons/' + run_name + '.txt' , np.transpose( [ angles_profile , datview ] ) , fmt='%.6f' ) 
	return [ angles_profile[3:] , datview [3:] , h[0] + h0 ]

def load_profile( file_name ):
	data = np.loadtxt( file_name )	
	return [ data[0,0] , data[0,1] , data[1,0] , data[3:,0] , data[3:,1] , data[3,0] , data[2,1] , data[2,0] , data[1,1] ]

def load_profile_border( file_name ):
	data = np.loadtxt( file_name )	
	return [ data[:,0] , data[:,1] ]

def onclick( event ):
	global newborderx
	global newbordery
	if( event.button == 3 ):
		newborderx = np.append( newborderx , event.xdata )
		newbordery = np.append( newbordery , event.ydata )
		plt.plot( newborderx , newbordery , 'r.' , markersize = 9 )
		plt.show()

def create_border_line( image , filename ):
	global newborderx
	global newbordery
	newborderx = np.zeros( 0 )
	newbordery = np.zeros( 0 )
	fig, ax = plt.subplots()
	ax.imshow( cv2.cvtColor( image , cv2.COLOR_BGR2RGB ) )
	messagebox.showinfo( title = None , message = "Instructions. Right-click to manage zoom. Left-click to add points." )
	cid = fig.canvas.mpl_connect( 'button_press_event' , onclick )
	plt.show()
	np.savetxt( 'ReferenceProfiles/' + filename + ".txt" , np.transpose( [ newborderx , newbordery ] ) , fmt = '%.4f' )
	cv2.imwrite( 'ReferenceProfiles/' + filename + ".png" , image )
	sort_index = np.argsort( newborderx )
	newborderx = newborderx[ sort_index ]
	newbordery = newbordery[ sort_index ]
	return [ newborderx , newbordery ]

def find_profile( compname , bordex , bordey , npositions , zoomfactor , angles_profiles , view_profiles , ninclination , inclination_1 , inclination_2 , orientation_1 , orientation_2 , stephor ):
	error_zoom = 1e10
	indexes_used = np.where( ( angles_profiles >= orientation_1 ) & ( angles_profiles <= orientation_2 ) )
	angles_profiles_used = angles_profiles[ indexes_used ]
	view_profiles_used = view_profiles[ indexes_used ]	
	[ list_datax , list_datay ] = normalize_inclined( bordex , -bordey , ninclination , inclination_1 , inclination_2 )
	inclination_steps = np.arange( inclination_1 , inclination_2 + ( inclination_2 - inclination_1  ) / ( ninclination - 1 ) / 2 , ( inclination_2 - inclination_1  ) / ( ninclination - 1 ) )
	candidates_x = np.zeros( 0 )
	candidates_y = np.zeros( 0 )
	candidates_z = np.zeros( 0 )
	candidates_sens = np.zeros( 0 )
	resolution = 1.0
	for i in range( len( list_datax ) ):
		datax = list_datax[ i ]
		datay = list_datay[ i ]
		stats_data = [ np.max( datay ) , np.polyfit( datax, datay, 1 )[0] , np.polyfit( datax, datay, 2 )[0] ]
		varver_border = np.max( datay ) - np.min( datay )
		varhor_border = np.max( datax ) - np.min( datax )
		asp_ratio_border = varver_border / varhor_border
		d1_border = np.polyfit( datax , datay , 1 )[0]
		d2_border = np.polyfit( datax , datay , 2 )[0]
		step_positions = ( np.max( angles_profiles_used ) - np.min( angles_profiles_used ) ) / npositions
		if( step_positions < 10 * stephor ):
			first_window = np.round( 10 * stephor / step_positions ) * step_positions
		else:
			first_window = step_positions
		view_init = np.arange( np.min( angles_profiles_used ) , np.max( angles_profiles_used ) , step_positions )
		windows = np.arange( first_window , np.max( angles_profiles_used ) - np.min( angles_profiles_used ) , step_positions )
		[ candidates_x , candidates_y , candidates_z , candidates_sens , resolution ] = candidates_initial( candidates_x , candidates_y , candidates_z , candidates_sens , inclination_steps[i] , stats_data , angles_profiles_used , view_profiles_used , view_init , windows , resolution , step_positions )
		print('First stage: Angle ' + str( np.round( inclination_steps[i] , 3 ) ) + ' finished.' )
	indexes_ok = np.where( candidates_sens <= resolution )
	candidates_x = candidates_x[ indexes_ok ]
	candidates_y = candidates_y[ indexes_ok ]
	candidates_z = candidates_z[ indexes_ok ]
	if( len( candidates_x ) > 0 ):
		[ candidates_x_zoom , candidates_y_zoom , candidates_inclination_zoom , error_zoom ] = candidates_zoom( candidates_x , candidates_y , candidates_z , step_positions , ( inclination_2 - inclination_1  ) / ( ninclination - 1 ) , zoomfactor , angles_profiles_used , view_profiles_used  , stats_data , resolution , bordex , bordey , error_zoom )
	if( compname != '' ):
		np.savetxt( 'CameraOrientations/' + compname + ".txt" , [ candidates_x_zoom[ -1 ] , candidates_y_zoom[ -1 ] , candidates_inclination_zoom[ -1 ] , npositions , zoomfactor , ninclination , inclination_1 , inclination_2 , orientation_1 , orientation_2 ] , fmt = '%.10f' )	
	return [ candidates_x_zoom[ -1 ] , candidates_y_zoom[ -1 ] , candidates_inclination_zoom[ -1 ] , error_zoom ]

def load_comparison( file_name ):
	data = np.loadtxt( file_name )
	return [ data[0] , data[1] , data[2] , data[3] , data[4] , data[5] , data[6] , data[7] , data[8] , data[9] ]

def renormalize( angles_profile , view_profile , init_view , range_view , inclination , borderx , bordery ):	
	cur_indexes = np.where( np.abs( angles_profile - ( init_view + range_view / 2 ) ) < range_view / 2 )
	datax = angles_profile[ cur_indexes ]
	datay = view_profile[ cur_indexes ]
	[ x0 , y0 ] = [ datax[0] , datay[0] ]
	new_xy = np.column_stack( ( datax - x0 , datay - y0 )) @ rotate( - inclination ) + [ x0 , y0 ]
	angles_current = ( new_xy[:,0] - np.min( new_xy[:,0] ) ) * ( np.max( borderx ) - np.min( borderx ) ) / ( np.max( new_xy[:,0] ) - np.min( new_xy[:,0] ) ) + np.min( borderx )
	view_current = ( new_xy[:,1] - np.min( new_xy[:,1] ) ) * ( np.max(borderx) - np.min(borderx) ) / ( np.max( new_xy[:,0] ) - np.min( new_xy[:,0] ) )
	f = interpolate.interp1d( angles_current , - view_current , fill_value = 'extrapolate' )
	datay_new = f( borderx )
	shift = np.mean( bordery ) - np.mean( datay_new )
	view_current = - view_current + shift
	return [ angles_current , view_current ]

def load_pix_height( foldername ):
	data = np.genfromtxt( foldername + '/Data.txt' , dtype = 'str' )
	for i in range( 0 , len( data ) - 1 ):
		cur_matrix = np.loadtxt( foldername + '/' + data[i,0] )
		if( i == 0 ):
			[ dimy, dimx ] = np.shape( cur_matrix )
			height_pix = np.zeros([ dimy , dimx , len( data ) ])
		height_pix[:,:,i] = cur_matrix
	list_planes = data[:,1]
	if( len( list_planes ) - 1 > 20 ):
		height_pix_stats = np.zeros([ dimy , dimx , 7 ])
		percentiles = [ 5 , 25 , 50 , 75 , 95 ]
	elif( len( list_planes ) - 1 > 10 ):
		height_pix_stats = np.zeros([ dimy , dimx , 5 ])
		percentiles = [ 25 , 50 , 75 ]		
	else:
		height_pix_stats = np.zeros([ dimy , dimx , 2 ])
		percentiles = []
	height_pix_stats[:,:,0] = np.loadtxt( foldername + '/Matrix_Average.txt' )
	height_pix_stats[:,:,1] = np.loadtxt( foldername + '/Matrix_Std.txt' )
	counter = 2
	for i in percentiles:
		height_pix_stats[:,:,counter] =  np.loadtxt( foldername + '/Matrix_Percentile_ ' + str( i ) + '.txt' ) 
		counter += 1
	height_pix_max = np.loadtxt( foldername + '/Matrix_Perpendicular.txt' )

	return [ height_pix , height_pix_stats , height_pix_max , list_planes ]

def angular_ranges( image , angles_profile , view_profile , init_view , range_view , inclination , borderx , bordery ):	
	f = interpolate.interp1d( angles_profile , view_profile )
	data_borders = [ f( init_view ) , f( init_view + range_view ) ]
	distang = np.sqrt( ( data_borders[1] - data_borders[0] ) * ( data_borders[1] - data_borders[0] ) + range_view * range_view )
	distpix = np.sqrt( ( borderx[-1] - borderx[0] ) * ( borderx[-1] - borderx[0] ) + ( bordery[-1] - bordery[0] ) * ( bordery[-1] - bordery[0] )  )
	dimx = image.shape[ 1 ]
	dimy = image.shape[ 0 ]
	height_pix_x = np.zeros([ dimy , dimx ])
	height_pix_y = np.zeros([ dimy , dimx ])
	for i in range( dimx ):
		for j in range( dimy ):
			[ height_pix_x[ j , i ] , height_pix_y[ j , i ] ] = [ i + 1 , j + 1 ]  @ rotate( - inclination )
	p1c = [ height_pix_x[ int( np.round( bordery[ 0 ] ) - 1 ), int( np.round( borderx[ 0 ] ) - 1 )], height_pix_y[ int( np.round( bordery[ 0 ] ) - 1 ),int( np.round(borderx[0])-1 )]]
	height_pix_q1 = ( height_pix_x - p1c[0] ) * distang / distpix + init_view
	height_pix_q2 = - ( height_pix_y - p1c[1] ) * distang / distpix + data_borders[ 0 ]
	Vec_x = np.cos( ( 90.0 - height_pix_q1 ) * np.pi / 180 )* np.cos( ( height_pix_q2 ) * np.pi / 180 )
	Vec_y = np.sin( ( 90.0 - height_pix_q1 ) * np.pi / 180 )* np.cos( ( height_pix_q2 ) * np.pi / 180 )
	Vec_z = np.sin( ( height_pix_q2 ) * np.pi / 180 )
	points = [ [ 0 , int( dimy / 2 ) ] , [ int( dimx / 2 ) , int( dimy / 2 ) ] , [ dimx - 1 , int( dimy / 2 )] ]
	data_angles = []
	for i in points:
		cur_vector = np.array( [ Vec_x[ i[1] , i[0] ] , Vec_y[ i[1] , i[0] ] , Vec_z[ i[1] , i[0] ] ] )
		cur_vector = cur_vector / np.sum( cur_vector * cur_vector )
		if( cur_vector[ 0 ] == 0 ):
			data_angles.append( 0.0 )
		elif( cur_vector[ 0 ] > 0 ):
			data_angles.append( 90.0 - np.arctan( cur_vector[ 1 ] / cur_vector[ 0 ] ) * 180.0 / np.pi )
		elif( cur_vector[ 1 ] >= 0 ):
			data_angles.append( -90.0 - np.arctan( cur_vector[ 1 ] / cur_vector[ 0 ] ) * 180.0 / np.pi )
		else:
			data_angles.append( 270.0 - np.arctan( cur_vector[ 1 ] / cur_vector[ 0 ] ) * 180.0 / np.pi )		
		data_angles.append( cur_vector[ 2 ] * 180 / np.pi )
	return data_angles

def pix_height( foldername , plane_1 , plane_2 , nplanes , maxhor , maxhei , minang , typemap , corx , cory , vent_corx , vent_cory , camheight , image , angles_profile , view_profile , init_view , range_view , inclination , borderx , bordery ):	
	if not os.path.exists( "PixelHeightConversion/" + foldername ):
		os.mkdir( "PixelHeightConversion/" + foldername )		
	str_name_data = "PixelHeightConversion/" + foldername + '/Data.txt'
	fdata = open( str_name_data , 'w' )
	f = interpolate.interp1d( angles_profile , view_profile )
	data_borders = [ f( init_view ) , f( init_view + range_view ) ]
	distang = np.sqrt( ( data_borders[1] - data_borders[0] ) * ( data_borders[1] - data_borders[0] ) + range_view * range_view )
	distpix = np.sqrt( ( borderx[-1] - borderx[0] ) * ( borderx[-1] - borderx[0] ) + ( bordery[-1] - bordery[0] ) * ( bordery[-1] - bordery[0] )  )
	dimx = image.shape[ 1 ]
	dimy = image.shape[ 0 ]
	height_pix_x = np.zeros([ dimy , dimx ])
	height_pix_y = np.zeros([ dimy , dimx ])
	for i in range( dimx ):
		for j in range( dimy ):
			[ height_pix_x[ j , i ] , height_pix_y[ j , i ] ] = [ i + 1 , j + 1 ]  @ rotate( - inclination )
	p1c = [ height_pix_x[ int( np.round( bordery[ 0 ] ) - 1 ), int( np.round( borderx[ 0 ] ) - 1 )], height_pix_y[ int( np.round( bordery[ 0 ] ) - 1 ),int( np.round(borderx[0])-1 )]]
	height_pix_q1 = ( height_pix_x - p1c[0] ) * distang / distpix + init_view
	height_pix_q2 = - ( height_pix_y - p1c[1] ) * distang / distpix + data_borders[ 0 ]
	Vec_x = np.cos( ( 90.0 - height_pix_q1 ) * np.pi / 180 )* np.cos( ( height_pix_q2 ) * np.pi / 180 )
	Vec_y = np.sin( ( 90.0 - height_pix_q1 ) * np.pi / 180 )* np.cos( ( height_pix_q2 ) * np.pi / 180 )
	Vec_z = np.sin( ( height_pix_q2 ) * np.pi / 180 )
	if( typemap == 1 ):
		camutm = utm.from_latlon( cory , corx )
		ventutm = utm.from_latlon( vent_cory , vent_corx , camutm[ 2 ] )
	else:
		camutm = [ corx , cory ]
		ventutm = [ vent_corx , vent_cory ]
	if( camutm[ 0 ] == ventutm[ 0 ] ):
		camang = 0.0
	else:
		dx =  camutm[ 0 ] - ventutm[ 0 ]
		dy =  camutm[ 1 ] - ventutm[ 1 ]
		camang = 90.0 - np.arctan( dy / dx ) * 180.0 / np.pi
	if( nplanes > 1 ):
		step_wind = np.abs( plane_2 - plane_1 ) / ( nplanes - 1 )
		winds = np.arange( min( plane_1 , plane_2 ) , max( plane_1 , plane_2 ) + step_wind / 2 , step_wind )		
	else:
		winds = [ plane_1 ]
	height_pix = np.zeros([ dimy , dimx , len( winds ) ])
	distance_pix = np.zeros([ dimy , dimx , len( winds ) ])
	counter = 0
	list_planes = []
	for wind in winds:
		if( wind >= 180 ):
			angwind = wind - 180.0
		else:
			angwind = wind
		diffang = abs( camang - wind )
		if( diffang >= minang and diffang <= 180 - minang ):
			Vec_pw = np.array( [ np.sin( ( wind - 90.0 ) * np.pi / 180.0 ) , np.cos( ( wind - 90.0 ) * np.pi / 180.0 )  , 0 ] )
			for i in range( dimx ):
				for j in range( dimy ):
					Vec_c = np.array( [ Vec_x[j,i] , Vec_y[j,i] , Vec_z[j,i] ] )
					Vec_c = Vec_c / np.sum( Vec_c * Vec_c )
					t = ( Vec_pw[0] * ( ventutm[0] - camutm[0] ) + Vec_pw[1] * ( ventutm[1] -  camutm[1] ) ) / ( Vec_pw[0] * Vec_c[0] + Vec_pw[1] * Vec_c[1] + Vec_pw[2] * Vec_c[2] )			
					dhor = np.sqrt( ( ventutm[0] - camutm[0] - t * Vec_c[ 0 ] ) * ( ventutm[0] - camutm[0] - t * Vec_c[ 0 ] ) + ( ventutm[1] - camutm[1] - t * Vec_c[ 1 ] ) * ( ventutm[1] - camutm[1] - t * Vec_c[ 1 ] ) )
					if( t >= 0 and dhor < maxhor ):
						distance_pix[j,i,counter] = dhor
						if( camheight + t * Vec_c[ 2 ] < maxhei ):
						    height_pix[j,i,counter] = camheight + t * Vec_c[ 2 ]
						else:
						    height_pix[j,i,counter] = np.nan
					else:
						distance_pix[j,i,counter] = np.nan
						height_pix[j,i,counter] = np.nan	
			str_name = "PixelHeightConversion/" + foldername + '/Matrix_' + str( round( counter + 1 ) ) + '.txt'
			np.savetxt( str_name , height_pix[:,:,counter] , fmt='%.2f' , delimiter=' ')
			np.savetxt( fdata , np.array( [ 'Matrix_' + str( round( counter + 1 ) ) + '.txt \t' + str( wind ) ] ) , fmt = "%s" )
			list_planes.append( str( wind ) ) 
			counter = counter + 1
			print( 'Direction ' + str( wind ) + ' finished.' )
		else:
			print( 'Direction ' + str( wind ) + ' skipped.' )
	height_pix = height_pix[ : , : , 0 : counter ]		
	height_pix_stats = np.nan
	if( counter > 1 ):
		if( counter > 20 ):
			height_pix_stats = np.zeros([ dimy , dimx , 7 ])
		elif( counter > 10 ):
			height_pix_stats = np.zeros([ dimy , dimx , 5 ])		
		else:
			height_pix_stats = np.zeros([ dimy , dimx , 2 ])
		height_pix_stats[ : , : , 0 ] = np.nanmean( height_pix , axis = 2 )
		str_name = "PixelHeightConversion/" + foldername + '/Matrix_Average.txt'
		np.savetxt( str_name , height_pix_stats[ : , : , 0 ] , fmt='%.2f', delimiter=' ')
		print( 'Average matrix finished.' )
		height_pix_stats[ : , : , 1 ] = np.nanstd( height_pix , axis = 2 )
		str_name = "PixelHeightConversion/" + foldername + '/Matrix_Std.txt'
		np.savetxt( str_name , height_pix_stats[ : , : , 1 ] , fmt='%.2f', delimiter=' ')
		print( 'Std matrix finished.' )
		if( counter > 10 ):
			percentiles = [ 25 , 50 , 75 ]
			if( counter > 20 ):
				percentiles = [ 5 , 25 , 50 , 75 , 95 ]
			counter2 = 2
			for i in percentiles:
				height_pix_stats[ : , : , counter2 ] = np.nanpercentile( height_pix , i , axis = 2 )
				str_name = "PixelHeightConversion/" + foldername + '/Matrix_Percentile_ ' + str( i ) + '.txt'
				np.savetxt( str_name , height_pix_stats[ : , : , counter2 ] , fmt='%.2f', delimiter=' ')
				print( 'Percentile ' + str( i ) + ' matrix finished.' )
				counter2 += 1
	if( camang >= 90 ):
		wind = camang - 90.0
	else:
		wind = camang + 90.0
	height_pix_max = np.zeros([ dimy , dimx ])
	Vec_pw = np.array( [ np.sin( ( wind - 90.0 ) * np.pi / 180.0 ) , np.cos( ( wind - 90.0 ) * np.pi / 180.0 )  , 0 ] )
	for i in range( dimx ):
		for j in range( dimy ):
			Vec_c = np.array( [ Vec_x[j,i] , Vec_y[j,i] , Vec_z[j,i] ] )
			t = ( Vec_pw[0] * ( ventutm[0] - camutm[0] ) + Vec_pw[1] * ( ventutm[1] -  camutm[1] ) ) / ( Vec_pw[0] * Vec_c[0] + Vec_pw[1] * Vec_c[1] + Vec_pw[2] * Vec_c[2] )			
			dhor = np.sqrt( ( ventutm[0] - camutm[0] - t * Vec_c[ 0 ] ) * ( ventutm[0] - camutm[0] - t * Vec_c[ 0 ] ) + ( ventutm[1] - camutm[1] - t * Vec_c[ 1 ] ) * ( ventutm[1] - camutm[1] - t * Vec_c[ 1 ] ) )
			if( t >= 0 and dhor < maxhor ):
				height_pix_max[j,i] = camheight + t * Vec_c[ 2 ]
			else:
				height_pix_max[j,i] = np.nan	
	str_name = "PixelHeightConversion/" + foldername + '/Matrix_Perpendicular.txt'
	np.savetxt( str_name , height_pix_max , fmt='%.2f' , delimiter=' ')
	np.savetxt( fdata , np.array( [ 'Perpendicular \t' + str( wind ) ] ) , fmt = "%s" ) 
	fdata.close()
	print( 'Direction perpendicular (' + str( round( wind , 2 ) ) +') finished.' )
	list_planes.append( str( wind ) ) 
	return [ height_pix , height_pix_stats , height_pix_max , list_planes ]

def set_levels( range_val ):
	mh = 50000.0
	if( range_val > 3000 ):
		return np.arange( 0 , mh , 500 )
	elif( range_val > 1500 ):
		levs = np.arange( 0 , mh , 250 )
	elif( range_val > 300 ):
		levs = np.arange( 0 , mh , 50 )
	elif( range_val > 150 ):
		levs = np.arange( 0 , mh , 25 )
	elif( range_val > 30 ):
		levs = np.arange( 0 , mh , 5 )
	else:
		levs = np.arange( 0 , mh , 1 )

def plot_pix_height( image , height_pix , height_pix_stats , height_pix_max , list_planes ):
	colortype = 'k'
	plt.imshow( cv2.cvtColor( image , cv2.COLOR_BGR2RGB ) , alpha = 0.3 )
	levs = set_levels( np.nanmax( height_pix_max ) - np.nanmin( height_pix_max ) )
	CS = plt.contour( height_pix_max , levels = levs , colors = colortype )
	plt.clabel( CS , fontsize = 10 , colors = colortype )
	plt.title( 'Perpendicular (' + str( round( float( list_planes[-1] ) , 3 ) ) + ')' )
	plt.show()
	for i in range( len( list_planes ) - 1 ):
		plt.imshow( cv2.cvtColor( image , cv2.COLOR_BGR2RGB ) , alpha = 0.3 )
		levs = set_levels( np.nanmax( height_pix[:,:,i] ) - np.nanmin( height_pix[:,:,i] ) )
		CS = plt.contour( height_pix[:,:,i] , levels = levs , colors = colortype )
		plt.clabel( CS , fontsize = 10 , colors = colortype )
		plt.title( str( round( float( list_planes[i] ) , 3 ) ) )
		plt.show()
	if( len( list_planes ) - 1 > 1 ):
		plt.imshow( cv2.cvtColor( image , cv2.COLOR_BGR2RGB ) , alpha = 0.3 )
		levs = set_levels( np.nanmax( height_pix_stats[:,:,0] ) - np.nanmin( height_pix_stats[:,:,0] ) )
		CS = plt.contour( height_pix_stats[:,:,0] , levels = levs , colors = colortype )
		plt.clabel( CS , fontsize = 10 , colors = colortype )
		plt.title( 'Average' )
		plt.show()
		plt.imshow( cv2.cvtColor( image , cv2.COLOR_BGR2RGB ) , alpha = 0.3 )
		if( np.max( height_pix_stats[:,:,1] ) > 500 ):
			levs = [ 30 , 50 , 100 , 500 , 1000 , 5000 ]
		elif( np.max( height_pix_stats[:,:,1] ) > 50 ):
			levs = [ 3 , 5 , 10 , 50 , 100 , 500 ]
		elif( np.max( height_pix_stats[:,:,1] ) > 5 ):
			levs = [ 0.3 , 0.5 , 1 , 5 , 10 , 50 ]	
		CS = plt.contour( height_pix_stats[:,:,1] , levels = levs , colors = colortype )
		plt.clabel( CS , fontsize = 10 , colors = colortype )
		plt.title( 'Std' )
		plt.show()
	if( len( list_planes ) - 1 > 10 ):
		percentiles = [ 25 , 50 , 75 ]
		if( len( list_planes ) - 1 > 20 ):
			percentiles = [ 5 , 25 , 50 , 75 , 95 ]
		counter2 = 2
		for i in percentiles:
			plt.imshow( cv2.cvtColor( image , cv2.COLOR_BGR2RGB ) , alpha = 0.3 )
			levs = set_levels( np.nanmax( height_pix_stats[:,:,counter2] ) - np.nanmin( height_pix_stats[:,:,counter2] ) )
			CS = plt.contour( height_pix_stats[:,:,counter2] , levels = levs , colors = colortype )
			plt.clabel( CS , fontsize = 10 , colors = colortype )
			plt.title( 'Percentile ' + str(i) )
			plt.show()
			counter2 += 1
    
def distance_two_points( lat1 , lat2 , lon1 , lon2 ):
	R = 6373.0
	lat1 = radians( lat1 )
	lon1 = radians( lon1 )
	lat2 = radians( lat2 )
	lon2 = radians( lon2 )
	dlon = lon2 - lon1
	dlat = lat2 - lat1
	a = sin( dlat / 2.0 ) ** 2.0 + cos( lat1 ) * cos( lat2 ) * sin( dlon / 2.0 ) ** 2.0
	c = 2.0 * atan2( sqrt( a ) , sqrt( 1.0 - a ) )
	return ( R * c ) * 1000.0
