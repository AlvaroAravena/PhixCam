from tkinter import END, Tk , ttk , Label , Button , W , E , Entry , OptionMenu , StringVar , DoubleVar , IntVar , Canvas , Scrollbar , VERTICAL , Frame , messagebox , filedialog
from PhixCam_Functions import get_profile_deg , get_profile_utm , import_map , read_map_deg , read_map_utm , plot_deg , plot_utm , load_profile, load_profile_border, load_comparison , find_profile, renormalize, pix_height , plot_pix_height , create_border_line, load_pix_height, angular_ranges
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import cv2
import warnings

class MainFrame:

	def __init__( self , master ):

		self.master = master
		self.image = -1
		self.source_dem_choice = 1
		self.topo_available = 0
		self.horizon_available = 0
		self.image_available = 0
		self.referenceprofile_available = 0
		self.comparison_available = 0
		self.matrices_available = 0
		self.source_dem = -1

		master.title( "PhixCam" )
		tab_parent = ttk.Notebook( master )
		tab1 = ttk.Frame( tab_parent )
		tab2 = ttk.Frame( tab_parent )
		tab3 = ttk.Frame( tab_parent )
		tab4 = ttk.Frame( tab_parent )

		for i in range( 4 ):
			tab1.columnconfigure( i , weight = 1 )
			tab2.columnconfigure( i , weight = 1 )
			tab3.columnconfigure( i , weight = 1 )
			tab4.columnconfigure( i , weight = 1 )
		input_runname = StringVar( master )
		input_runname.set( "HorizonProfile" )
		input_toponame = StringVar( master )
		input_toponame.set( "Topography" )
		input_points = IntVar( master )
		input_points.set( 1000 )
		input_lon = DoubleVar( master )
		input_lon.set( -71.973858 )
		input_lat = DoubleVar( master )
		input_lat.set( -39.276109 )
		input_east = DoubleVar( master )
		input_east.set( 630000.00 )
		input_north = DoubleVar( master )
		input_north.set( 7415000.00 )
		input_ang1 = DoubleVar( master )
		input_ang1.set( 150.0 )
		input_ang2 = DoubleVar( master )
		input_ang2.set( 180.0 )
		input_h0 = DoubleVar( master )
		input_h0.set( 15.0 )
		input_d_step = DoubleVar( master )
		input_d_step.set( 100.0 )
		input_max_distance = DoubleVar( master )
		input_max_distance.set( 100.0 )
		var_dem = StringVar( master )
		var_dem.set( "STRM 30 m" )
		input_compname = StringVar( master )
		input_compname.set( "CameraOrientation_1" )
		input_npositions = IntVar( master )
		input_npositions.set( 100 )
		input_ninclination = IntVar( master )
		input_ninclination.set( 11 )
		input_zoomfactor = IntVar( master )
		input_zoomfactor.set( 10 )
		input_orientation_1 = DoubleVar( master )
		input_orientation_1.set( 0.0 )
		input_orientation_2 = DoubleVar( master )
		input_orientation_2.set( 360.0 )
		input_inclination_1 = DoubleVar( master )
		input_inclination_1.set( -5.0 )
		input_inclination_2 = DoubleVar( master )
		input_inclination_2.set( 5.0 )
		input_savepoints = StringVar( master )
		input_savepoints.set( "Camera_1" )
		input_vent_lon = DoubleVar( master )
		input_vent_lon.set( -71.93927 )
		input_vent_lat = DoubleVar( master )
		input_vent_lat.set( -39.42087 )
		input_vent_east = DoubleVar( master )
		input_vent_east.set( 635000.00 )
		input_vent_north = DoubleVar( master )
		input_vent_north.set( 7420000.00 )
		input_plane_1 = DoubleVar( master )
		input_plane_1.set( 0.0 )
		input_plane_2 = DoubleVar( master )
		input_plane_2.set( 180.0 )
		input_nplanes = DoubleVar( master )
		input_nplanes.set( 1 )
		input_save_folder = StringVar( master )
		input_save_folder.set( "Results_1" )
		input_maxhor = DoubleVar( master )
		input_maxhor.set( 10000.0 )
		input_maxhei = DoubleVar( master )
		input_maxhei.set( 30000.0 )
		input_minang = DoubleVar( master )
		input_minang.set( 30.0 )
																			
		self.label_toposource = Label( tab1 , text = "DEM source" )
		self.label_toposource.grid( row = 0 , column = 0 , columnspan = 1 , sticky = W )
		self.toposource_entry = OptionMenu( tab1 , var_dem , "STRM 30 m" , "Input DEM (utm)" , "Input DEM (latlon)" , command = self.opt_dem )
		self.toposource_entry.grid( row = 0 , column = 1 , columnspan = 3 , sticky = W + E )
		self.label_runname = Label( tab1 , text = "Horizon Profile Name (when saving)" )
		self.label_runname.grid( row = 1 , column = 0 , columnspan = 1 , sticky = W )
		self.runname_entry = Entry( tab1 , textvariable = input_runname )
		self.runname_entry.grid( row = 1 , column = 1 , columnspan = 1 , sticky = E )
		self.label_toponame = Label( tab1 , text = "Topography Name (when saving)" )
		self.label_toponame.grid( row = 1 , column = 2 , columnspan = 1 , sticky = W )
		self.toponame_entry = Entry( tab1 , textvariable = input_toponame )
		self.toponame_entry.grid( row = 1 , column = 3 , columnspan = 1 , sticky = E )
		self.label_points = Label( tab1 , text = "Number of points" )
		self.label_points.grid( row = 2 , column = 0 , columnspan = 1 , sticky = W )
		self.points_entry = Entry( tab1 , textvariable = input_points )
		self.points_entry.grid( row = 2 , column = 1 , columnspan = 1 , sticky = E )
		self.label_h0 = Label( tab1 , text = "Camera Height (m)" )
		self.label_h0.grid( row = 2 , column = 2 , columnspan = 1 , sticky = W )
		self.h0_entry = Entry( tab1 , textvariable = input_h0 )
		self.h0_entry.grid( row = 2 , column = 3 , columnspan = 1 , sticky = E )
		self.label_lon = Label( tab1 , text = "Camera Longitude (deg)" )
		self.label_lon.grid( row = 3 , column = 0 , columnspan = 1 , sticky = W )
		self.lon_entry = Entry( tab1 , textvariable = input_lon )
		self.lon_entry.grid( row = 3 , column = 1 , columnspan = 1 , sticky = E )
		self.label_lat = Label( tab1 , text = "Camera Latitude (deg)" )
		self.label_lat.grid( row = 3 , column = 2 , columnspan = 1 , sticky = W )
		self.lat_entry = Entry( tab1 , textvariable = input_lat )
		self.lat_entry.grid( row = 3  , column = 3 , columnspan = 1 , sticky = E )
		self.label_east = Label( tab1 , text = "Camera East Position (m)" )
		self.label_east.grid( row = 4 , column = 0 , columnspan = 1 , sticky = W )
		self.east_entry = Entry( tab1 , textvariable = input_east )
		self.east_entry.grid( row = 4 , column = 1 , columnspan = 1 , sticky = E )
		self.label_north = Label( tab1 , text = "Camera North Position (m)" )
		self.label_north.grid( row = 4 , column = 2 , columnspan = 1 , sticky = W )
		self.north_entry = Entry( tab1 , textvariable = input_north )
		self.north_entry.grid( row = 4  , column = 3 , columnspan = 1 , sticky = E )
		self.label_ang1 = Label( tab1 , text = "Initial Direction (deg)" )
		self.label_ang1.grid( row = 5 , column = 0 , columnspan = 1 , sticky = W )
		self.ang1_entry = Entry( tab1 , textvariable = input_ang1 )
		self.ang1_entry.grid( row = 5 , column = 1 , columnspan = 1 , sticky = E )
		self.label_ang2 = Label( tab1 , text = "Final Direction (deg)" )
		self.label_ang2.grid( row = 5 , column = 2 , columnspan = 1 , sticky = W )
		self.ang2_entry = Entry( tab1 , textvariable = input_ang2 )
		self.ang2_entry.grid( row = 5 , column = 3 , columnspan = 1 , sticky = E )
		self.label_d = Label( tab1 , text = "Distance Step (m)" )
		self.label_d.grid( row = 6 , column = 0 , columnspan = 1 , sticky = W )
		self.d_step_entry = Entry( tab1 , textvariable = input_d_step )
		self.d_step_entry.grid( row = 6 , column = 1 , columnspan = 1 , sticky = E )
		self.label_maxdist = Label( tab1 , text = "Maximum distance (km)" )
		self.label_maxdist.grid( row = 6 , column = 2 , columnspan = 1 , sticky = W )
		self.max_distance_entry = Entry( tab1 , textvariable = input_max_distance )
		self.max_distance_entry.grid( row = 6 , column = 3 , columnspan = 1 , sticky = E )
		self.bot_load_topo = Button( tab1 , text = "Load Topography" , command = self.load_topo )
		self.bot_load_topo.grid( row = 7, column = 0 , columnspan = 2 , sticky = W + E)
		self.bot_plot_topo = Button( tab1 , text = "Plot Topography" , command = self.plot_topo )
		self.bot_plot_topo.grid( row = 7 , column = 2 , columnspan = 2 , sticky = W + E)
		self.bot_create_toptrace = Button( tab1 , text = "Create Horizon Profile" , command = self.create_toptrace )
		self.bot_create_toptrace.grid( row = 8, column = 0 , columnspan = 2 , sticky = W + E)
		self.bot_load_toptrace = Button( tab1 , text = "Load Horizon Profile" , command = self.load_toptrace )
		self.bot_load_toptrace.grid( row = 8 , column = 2 , columnspan = 2 , sticky = W + E)
		self.bot_plot_toptrace = Button( tab1 , text = "Plot Horizon Profile" , command = self.plot_toptrace )
		self.bot_plot_toptrace.grid( row = 9 , column = 0 , columnspan = 4 , sticky = W + E)

		self.label_refim = Label( tab2 , text = "Reference Image" )
		self.label_refim.grid( row = 0 , column = 0 , columnspan = 4 , sticky = W + E )
		self.bot_load_image = Button( tab2 , text = "Load Image" , command = self.load_image )
		self.bot_load_image.grid( row = 1 , column = 0 , columnspan = 2 , sticky = W + E )
		self.bot_show_image = Button( tab2 , text = "Show Image" , command = self.show_image )
		self.bot_show_image.grid( row = 1 , column = 2 , columnspan = 2 , sticky = W + E )
		self.label_refimpro = Label( tab2 , text = "Reference Profile" )
		self.label_refimpro.grid( row = 2 , column = 0 , columnspan = 4 , sticky = W + E )
		self.label_savepoints = Label( tab2 , text = "Reference Profile Name (when saving)" )
		self.label_savepoints.grid( row = 3 , column = 0 , columnspan = 2 , sticky = W )
		self.savepoints_entry = Entry( tab2 , textvariable =  input_savepoints )
		self.savepoints_entry.grid( row = 3 , column = 2 , columnspan = 2 , sticky = W + E )
		self.bot_create_border = Button( tab2 , text = "Create Profile" , command = self.create_border )
		self.bot_create_border.grid( row = 4 , column = 0 , columnspan = 2 , sticky = W + E )
		self.bot_load_border = Button( tab2 , text = "Load Profile" , command = self.load_border )
		self.bot_load_border.grid( row = 4 , column = 2 , columnspan = 2 , sticky = W + E )
		self.bot_plot_border = Button( tab2 , text = "Plot Profile" , command = self.plot_border )
		self.bot_plot_border.grid( row = 5 , column = 0 , columnspan = 4 , sticky = W + E )

		self.label_compname = Label( tab3 , text = "Comparison Name (when saving)" )
		self.label_compname.grid( row = 0 , column = 0 , columnspan = 2 , sticky = W )
		self.compname_entry = Entry( tab3 , textvariable = input_compname )
		self.compname_entry.grid( row = 0 , column = 2 , columnspan = 2 , sticky = E )
		self.label_npositions = Label( tab3 , text = "Number of Tested Initial Positions" )
		self.label_npositions.grid( row = 1 , column = 0 , columnspan = 2 , sticky = W )
		self.npositions_entry = Entry( tab3 , textvariable = input_npositions )
		self.npositions_entry.grid( row = 1 , column = 2 , columnspan = 2 , sticky = E )
		self.label_ninclination = Label( tab3 , text = "Number of Tested Camera Rotation Angles (roll)" )
		self.label_ninclination.grid( row = 2 , column = 0 , columnspan = 2 , sticky = W )
		self.ninclination_entry = Entry( tab3 , textvariable = input_ninclination )
		self.ninclination_entry.grid( row = 2 , column = 2 , columnspan = 2 , sticky = E )
		self.label_zoomfactor = Label( tab3 , text = "Zoom Factor" )
		self.label_zoomfactor.grid( row = 3 , column = 0 , columnspan = 2 , sticky = W )
		self.zoomfactor_entry = Entry( tab3 , textvariable = input_zoomfactor )
		self.zoomfactor_entry.grid( row = 3 , column = 2 , columnspan = 2 , sticky = E )
		self.label_orientation_1 = Label( tab3 , text = "Initial Orientation (deg)" )
		self.label_orientation_1.grid( row = 4 , column = 0 , columnspan = 2 , sticky = W )
		self.orientation_1_entry = Entry( tab3 , textvariable = input_orientation_1 )
		self.orientation_1_entry.grid( row = 4 , column = 2 , columnspan = 2 , sticky = E )
		self.label_orientation_2 = Label( tab3 , text = "Final Orientation (deg)" )
		self.label_orientation_2.grid( row = 5 , column = 0 , columnspan = 2 , sticky = W )
		self.orientation_2_entry = Entry( tab3 , textvariable = input_orientation_2 )
		self.orientation_2_entry.grid( row = 5 , column = 2 , columnspan = 2 , sticky = E )	
		self.label_inclination_1 = Label( tab3 , text = "Initial Rotation Angle (roll; deg)" )
		self.label_inclination_1.grid( row = 6 , column = 0 , columnspan = 2 , sticky = W )
		self.inclination_1_entry = Entry( tab3 , textvariable = input_inclination_1 )
		self.inclination_1_entry.grid( row = 6 , column = 2 , columnspan = 2 , sticky = E )
		self.label_inclination_2 = Label( tab3 , text = "Final Rotation Angle (roll; deg)" )
		self.label_inclination_2.grid( row = 7 , column = 0 , columnspan = 2 , sticky = W )
		self.inclination_2_entry = Entry( tab3 , textvariable = input_inclination_2 )
		self.inclination_2_entry.grid( row = 7 , column = 2 , columnspan = 2 , sticky = E )
		self.bot_load_comp_profile = Button( tab3 , text = "Load Comparison" , command = self.load_comparison )
		self.bot_load_comp_profile.grid( row = 8 , column = 0 , columnspan = 2 , sticky = W + E )
		self.bot_comp_profile = Button( tab3 , text = "Compare Profiles" , command = self.compare_profiles )
		self.bot_comp_profile.grid( row = 8 , column = 2 , columnspan = 2 , sticky = W + E )
		self.bot_plot_comp_profile = Button( tab3 , text = "Plot Comparison" , command = self.plot_compare_profiles )
		self.bot_plot_comp_profile.grid( row = 9 , column = 0 , columnspan = 2 , sticky = W + E )
		self.bot_improve_profile = Button( tab3 , text = "Set Inputs to Improve Comparison" , command = self.improve_profiles )
		self.bot_improve_profile.grid( row = 9 , column = 2 , columnspan = 2 , sticky = W + E )

		self.label_save_folder = Label( tab4 , text = "Save Folder Name" )
		self.label_save_folder.grid( row = 0 , column = 0 , columnspan = 2 , sticky = W )
		self.save_folder_entry = Entry( tab4 , textvariable = input_save_folder )
		self.save_folder_entry.grid( row = 0 , column = 2 , columnspan = 2 , sticky = E )
		self.label_vent_lon = Label( tab4 , text = "Vent Longitude (deg)" )
		self.label_vent_lon.grid( row = 1 , column = 0 , columnspan = 2 , sticky = W )
		self.vent_lon_entry = Entry( tab4 , textvariable = input_vent_lon )
		self.vent_lon_entry.grid( row = 1 , column = 2 , columnspan = 2 , sticky = E )
		self.label_vent_lat = Label( tab4 , text = "Vent Latitude (deg)" )
		self.label_vent_lat.grid( row = 2 , column = 0 , columnspan = 2 , sticky = W )
		self.vent_lat_entry = Entry( tab4 , textvariable = input_vent_lat )
		self.vent_lat_entry.grid( row = 2  , column = 2 , columnspan = 2 , sticky = E )
		self.label_vent_east = Label( tab4 , text = "Vent East Position (m)" )
		self.label_vent_east.grid( row = 3 , column = 0 , columnspan = 2 , sticky = W )
		self.vent_east_entry = Entry( tab4 , textvariable = input_vent_east )
		self.vent_east_entry.grid( row = 3 , column = 2 , columnspan = 2 , sticky = E )
		self.label_vent_north = Label( tab4 , text = "Vent North Position (m)" )
		self.label_vent_north.grid( row = 4 , column = 0 , columnspan = 2 , sticky = W )
		self.vent_north_entry = Entry( tab4 , textvariable = input_vent_north )
		self.vent_north_entry.grid( row = 4  , column = 2 , columnspan = 2 , sticky = E )
		self.label_plane_1 = Label( tab4 , text = "Initial Plane Orientation (deg)" )
		self.label_plane_1.grid( row = 5 , column = 0 , columnspan = 2 , sticky = W )
		self.plane_1_entry = Entry( tab4 , textvariable = input_plane_1 )
		self.plane_1_entry.grid( row = 5  , column = 2 , columnspan = 2 , sticky = E )
		self.label_plane_2 = Label( tab4 , text = "Final Plane Orientation (deg)" )
		self.label_plane_2.grid( row = 6 , column = 0 , columnspan = 2 , sticky = W )
		self.plane_2_entry = Entry( tab4 , textvariable = input_plane_2 )
		self.plane_2_entry.grid( row = 6  , column = 2 , columnspan = 2 , sticky = E )
		self.label_nplanes = Label( tab4 , text = "Number of Steps" )
		self.label_nplanes.grid( row = 7 , column = 0 , columnspan = 2 , sticky = W )
		self.nplanes_entry = Entry( tab4 , textvariable = input_nplanes )
		self.nplanes_entry.grid( row = 7  , column = 2 , columnspan = 2 , sticky = E )
		self.label_maxhor = Label( tab4 , text = "Max. horizontal distance from vent (m)" )
		self.label_maxhor.grid( row = 8 , column = 0 , columnspan = 2 , sticky = W )
		self.maxhor_entry = Entry( tab4 , textvariable = input_maxhor )
		self.maxhor_entry.grid( row = 8  , column = 2 , columnspan = 2 , sticky = E )
		self.label_maxhei = Label( tab4 , text = "Max. vertical height (m a.s.l.)" )
		self.label_maxhei.grid( row = 9 , column = 0 , columnspan = 2 , sticky = W )
		self.maxhei_entry = Entry( tab4 , textvariable = input_maxhei )
		self.maxhei_entry.grid( row = 9  , column = 2 , columnspan = 2 , sticky = E )
		self.label_minang = Label( tab4 , text = "Min. angular diff. between plane and camera orientation (deg)" )
		self.label_minang.grid( row = 10 , column = 0 , columnspan = 2 , sticky = W )		
		self.minang_entry = Entry( tab4 , textvariable = input_minang )
		self.minang_entry.grid( row = 10  , column = 2 , columnspan = 2 , sticky = E )
		self.bot_loadpixheight = Button( tab4 , text = "Load matrices" , command = self.load_pixheight )
		self.bot_loadpixheight.grid( row = 11 , column = 0 , columnspan = 1 , sticky = W + E)
		self.bot_pixheight = Button( tab4 , text = "Create matrices" , command = self.create_pixheight )
		self.bot_pixheight.grid( row = 11 , column = 1 , columnspan = 1 , sticky = W + E)
		self.bot_plot_pixheight = Button( tab4 , text = "Plot matrices" , command = self.plot_pixheight )
		self.bot_plot_pixheight.grid( row = 11 , column = 2 , columnspan = 2 , sticky = W + E)
		
		tab_parent.add( tab1 , text = "Horizon Profile" )
		tab_parent.add( tab2 , text = "Reference Profile" )
		tab_parent.add( tab3 , text = "Georeference Image" )
		tab_parent.add( tab4 , text = "Pixel Height Conversion" )
		tab_parent.pack( expand = 1 , fill = 'both' )

		self.enabled_disabled()

	def opt_dem( self , opt ):
		if( opt == "STRM 30 m" ):
			self.source_dem_choice = 1
		elif( opt == "Input DEM (utm)" ):
			self.source_dem_choice = 2
		else:
			self.source_dem_choice = 3
		self.enabled_disabled()

	def show_image( self ):
		plt.figure( 1 )
		plt.imshow( cv2.cvtColor( self.image , cv2.COLOR_BGR2RGB ) )
		plt.show()

	def load_topo( self ):
		prev_source_dem = self.source_dem
		self.source_dem = self.source_dem_choice
		if( self.source_dem == 1 ):
			self.toponame = self.toponame_entry.get()
			loadlon = float( self.lon_entry.get() )
			loadlat = float( self.lat_entry.get() )
			[ self.lon1 , self.lon2 , self.lat1 , self.lat2 , self.Cities , self.Topography , self.Topography_Sea , self.cells_lon , self.cells_lat ] = import_map( self.toponame , loadlon , loadlat )
		elif( self.source_dem == 2 ):
			file_path = filedialog.askopenfilename( initialdir = "Topographies" , title = "Select image file", filetypes = [("Text files", "*.asc")] )
			[ self.Topography , self.Topography_Sea , self.n_north , self.n_east , self.cellsize , self.east_cor , self.north_cor ] = read_map_utm( 'Topographies' , file_path )
		elif( self.source_dem == 3 ):
			file_path = filedialog.askopenfilename( initialdir = "Topographies" , title = "Select image file", filetypes = [("Text files", "*.txt")] )
			[ self.lon1 , self.lon2 , self.lat1 , self.lat2 , self.Cities , self.Topography , self.Topography_Sea , self.cells_lon , self.cells_lat ] = read_map_deg( 'Topographies' , file_path )
		if( prev_source_dem > -1 ):
			if( ( self.source_dem in [ 1 , 3 ] and prev_source_dem == 2 ) or ( self.source_dem == 2 and prev_source_dem in [ 1 , 3 ] ) ):
				self.horizon_available = 0
				messagebox.showinfo( title = None , message = "Loaded topography DEM type is not consistent with horizon profile. Imported horizon profile was discarded." )
		self.topo_available = 1
		self.enabled_disabled()
		messagebox.showinfo( title = None , message = "Topography loaded successfully" )

	def plot_topo( self ):
		plotang1 = float( self.ang1_entry.get() )
		plotang2 = float( self.ang2_entry.get() )
		if( self.source_dem in [ 1 , 3 ] ):
			plotlon = float( self.lon_entry.get() )
			plotlat = float( self.lat_entry.get() )
			plot_deg( plotlon , plotlat , self.lon1 , self.lon2 , self.lat1 , self.lat2 , self.Cities , self.Topography , self.Topography_Sea , self.cells_lon , self.cells_lat , plotang1 , plotang2 )
		else:
			ploteast = float( self.east_entry.get() )
			plotnorth = float( self.north_entry.get() )
			plot_utm( ploteast , plotnorth , self.east_cor , self.north_cor , self.n_east , self.n_north , self.cellsize , self.Topography , self.Topography_Sea , plotang1 , plotang2 )

	def create_toptrace( self ):
		self.runname = self.runname_entry.get()
		numpoints = float( self.points_entry.get() )
		self.ang1 = float( self.ang1_entry.get() )
		self.ang2 = float( self.ang2_entry.get() )
		self.h0 = float( self.h0_entry.get() )
		self.d_step = float( self.d_step_entry.get() )
		self.max_distance = float( self.max_distance_entry.get() )
		if( self.source_dem in [ 1 , 3 ] ):
			self.lon = float( self.lon_entry.get() )
			self.lat = float( self.lat_entry.get() )
			[ self.angles_profile , self.view_profile , self.camheight ] = get_profile_deg( self.runname , numpoints , self.lon , self.lat ,  self.lon1 , self.lon2 , self.lat1 , self.lat2 , self.Cities , self.Topography , self.Topography_Sea , self.cells_lon , self.cells_lat , self.ang1 , self.ang2 , self.h0 , self.d_step , self.max_distance )
		else:
			self.east = float( self.east_entry.get() )
			self.north = float( self.north_entry.get() )
			[ self.angles_profile , self.view_profile , self.camheight ] = get_profile_utm( self.runname , numpoints , self.east , self.north , self.east_cor , self.north_cor , self.cellsize , self.Topography , self.n_east , self.n_north , self.ang1 , self.ang2 , self.h0 , self.d_step , self.max_distance )
		self.stephor = self.angles_profile[ 1 ] - self.angles_profile[ 0 ]
		self.horizon_available = 1
		self.enabled_disabled()
		messagebox.showinfo( title = None , message = "Profile saved successfully" )

	def load_toptrace( self ):
		file_path = filedialog.askopenfilename( initialdir = "Horizons" , title = "Select topography profile file" , filetypes=[("Text files", "*.txt")])
		prev_source_dem = self.source_dem
		[ dn_lon , dn_lat , self.camheight , self.angles_profile , self.view_profile , self.ang1 , self.ang2 , self.h0 , dn_source_dem ] = load_profile( file_path )
		if( dn_source_dem == 1.0 ):
			[ self.lon , self.lat ] = [ dn_lon , dn_lat ]
			self.source_dem = 1
			self.source_dem_choice = 1
		else:
			[ self.east , self.north ] = [ dn_lon , dn_lat ]
			self.source_dem = 2
			self.source_dem_choice = 2
		self.enabled_disabled()
		if( prev_source_dem > -1 ):
			if( ( self.source_dem in [ 1 , 3 ] and prev_source_dem == 2 ) or ( self.source_dem == 2 and prev_source_dem in [ 1 , 3 ] ) ):
				self.topo_available = 0
				messagebox.showinfo( title = None , message = "Loaded topography DEM type is not consistent with horizon profile. Imported topography was discarded." )
		self.stephor = self.angles_profile[ 1 ] - self.angles_profile[ 0 ]
		self.horizon_available = 1
		if( self.source_dem in [ 1 , 3 ] ):
			self.lon_entry.delete( 0 , len( self.lon_entry.get() ) )
			self.lon_entry.insert( 0 , str( self.lon ) )
			self.lat_entry.delete( 0 , len( self.lat_entry.get() ) )
			self.lat_entry.insert( 0 , str( self.lat ) )
		else:
			self.east_entry.delete( 0 , len( self.east_entry.get() ) )
			self.east_entry.insert( 0 , str( self.east ) )
			self.north_entry.delete( 0 , len( self.north_entry.get() ) )
			self.north_entry.insert( 0 , str( self.north ) )
		self.h0_entry.delete( 0 , len( self.h0_entry.get() ) )
		self.h0_entry.insert( 0 , str( self.h0 ) )
		self.ang1_entry.delete( 0 , len( self.ang1_entry.get() ) )
		self.ang1_entry.insert( 0 , str( self.ang1 ) )
		self.ang2_entry.delete( 0 , len( self.ang2_entry.get() ) )
		self.ang2_entry.insert( 0 , str( self.ang2 ) )
		self.enabled_disabled()
		messagebox.showinfo( title = None , message = "Profile loaded successfully" )

	def plot_toptrace( self ):
		plt.figure( 2 )
		plt.plot( self.angles_profile , self.view_profile , "k-" )
		ax = plt.gca()
		ax.set_aspect('equal', adjustable='box')
		ax.set_xlabel('Horizontal angle (deg)')
		ax.set_ylabel('Vertical angle (deg)')
		plt.show()
		
	def load_image( self ):
		file_path = filedialog.askopenfilename( title="Select image file", filetypes=[("Image files", "*.png *.jpg ")])
		self.image = cv2.imread( file_path )
		self.image_available = 1
		self.enabled_disabled()
		messagebox.showinfo( title = None , message = "Reference image loaded successfully" )
		
	def create_border( self ):
		[ self.borderx , self.bordery ] = create_border_line( self.image , self.savepoints_entry.get() )
		self.referenceprofile_available = 1
		self.enabled_disabled()
		messagebox.showinfo( title = None , message = "Reference profile created successfully" )
						
	def load_border( self ):
		file_path = filedialog.askopenfilename( initialdir = "ReferenceProfiles" , title = "Select reference profile file", filetypes=[("Text files", "*.txt")])
		[ self.borderx , self.bordery ] = load_profile_border( file_path )
		self.image =cv2.imread( file_path.replace( 'txt' , 'png' ) )
		self.referenceprofile_available = 1
		self.image_available = 1
		self.enabled_disabled()
		messagebox.showinfo( title = None , message = "Reference profile loaded successfully" )
		
	def plot_border( self ):
		plt.figure( 3 )
		plt.imshow( cv2.cvtColor( self.image , cv2.COLOR_BGR2RGB ) )
		plt.plot( self.borderx , self.bordery , 'r.' , markersize = 9 )
		plt.show()
	
	def load_comparison( self ):
		file_path = filedialog.askopenfilename( initialdir = "CameraOrientations" , title="Select comparison file", filetypes=[("Text files", "*.txt")])
		[ self.fix_ang1 , self.fix_ang2 , self.fix_inclination , self.npositions , self.zoomfactor , self.ninclination , self.inclination_1 , self.inclination_2 , self.orientation_1 , self.orientation_2 ] = load_comparison( file_path )
		self.comparison_available = 1
		self.enabled_disabled()
		messagebox.showinfo( title = None , message = "Georeferenciation loaded successfully" )

	def compare_profiles( self ):
		self.npositions = float( self.npositions_entry.get() )
		self.ninclination = float( self.ninclination_entry.get() )
		self.zoomfactor = float( self.zoomfactor_entry.get() )
		self.orientation_1 = float( self.orientation_1_entry.get() )
		self.orientation_2 = float( self.orientation_2_entry.get() )
		self.inclination_1 = float( self.inclination_1_entry.get() )
		self.inclination_2 = float( self.inclination_2_entry.get() )
		self.compname = self.compname_entry.get()
		[ self.fix_ang1 , self.fix_ang2 , self.fix_inclination , fit_value ] = find_profile( self.compname , self.borderx , self.bordery , self.npositions , self.zoomfactor , self.angles_profile , self.view_profile , self.ninclination , self.inclination_1 , self.inclination_2 , self.orientation_1 , self.orientation_2 , self.stephor )
		data_angular = angular_ranges( self.image , self.angles_profile , self.view_profile , self.fix_ang1 , self.fix_ang2 , self.fix_inclination , self.borderx , self.bordery )
		print(' ')
		if( np.sqrt(fit_value / len( self.borderx ) ) > 0.01 ):
		    print('Error measure is too high, please check the range used for roll rotation and the number of tested values.')
		    print(' ')
		print( 'Best fit' )
		print( 'Angle of first point: ' + str( round( self.fix_ang1 , 2 ) ) + ' deg.' )
		print( 'Angle of last point: ' + str( round( self.fix_ang1 + self.fix_ang2 , 2 ) ) + ' deg.' )
		print( 'Angular distance between first and last point: ' + str( round( self.fix_ang2 , 2 ) )  + ' deg.' )
		print( 'Roll Rotation: ' + str( round( self.fix_inclination , 2 ) ) + ' deg.' )
		print( 'Tilt Rotation at Image Center: ' + str( round( data_angular[ 3 ] , 2 ) )  + ' deg.' )
		print( 'Angle of left border: ' + str( round( data_angular[ 0 ] , 2 ) )  + ' deg.' )
		print( 'Angle of right border: ' + str( round( data_angular[ 4 ] , 2 ) )  + ' deg.' )
		print( 'Error measure: ' + str( round( np.sqrt(fit_value / len( self.borderx ) ) , 5 ) )  )
		self.comparison_available = 1
		self.enabled_disabled()
		messagebox.showinfo( title = None , message = "Georeferenciation performed successfully" )

	def improve_profiles( self ):
		step_positions = ( self.orientation_2 - self.orientation_1 ) / ( self.npositions - 1 )
		self.orientation_1_entry.delete( 0 , len( self.orientation_1_entry.get() ) )
		self.orientation_1_entry.insert( 0 , str( round( self.fix_ang1 - 5 * step_positions , 3 ) ) )		
		self.orientation_2_entry.delete( 0 , len( self.orientation_2_entry.get() ) )
		self.orientation_2_entry.insert( 0 , str( round( self.fix_ang1 + self.fix_ang2 + 5 * step_positions , 3 ) ) )
		step_inclinations = ( self.inclination_2 - self.inclination_1 ) / ( self.ninclination - 1 )	
		self.inclination_1_entry.delete( 0 , len( self.inclination_1_entry.get() ) )
		self.inclination_1_entry.insert( 0 , str( round( self.fix_inclination - step_inclinations , 3 ) ))	
		self.inclination_2_entry.delete( 0 , len( self.inclination_2_entry.get() ) )
		self.inclination_2_entry.insert( 0 , str( round( self.fix_inclination + step_inclinations , 3 ) ) )		
						
	def plot_compare_profiles( self ):
		plt.figure( 4 )
		plt.imshow( cv2.cvtColor( self.image , cv2.COLOR_BGR2RGB ) )
		plt.plot( self.borderx , self.bordery , 'r.' , markersize = 9 )
		[ nor_angles_profile , nor_view_profile ] = renormalize( self.angles_profile , self.view_profile , self.fix_ang1 , self.fix_ang2 , self.fix_inclination , self.borderx , self.bordery )
		ax = plt.gca()
		ax.set_aspect( 'equal', adjustable = 'box' )
		plt.plot( nor_angles_profile , nor_view_profile , "b:" , linewidth = 2 )
		plt.show()

	def load_pixheight( self ):
		folder_path = filedialog.askdirectory( initialdir = "PixelHeightConversion" , title = "Select folder of conversion matrices" )
		[ self.pix_matrices , self.pix_matrices_stats , self.pix_matrices_max , self.list_planes ] = load_pix_height( folder_path )
		self.matrices_available = 1
		self.enabled_disabled()
		messagebox.showinfo( title = None , message = "Matrices loaded successfully" )
		
	def create_pixheight( self ):
		self.foldername = self.save_folder_entry.get()
		self.plane_1 = float( self.plane_1_entry.get() )
		self.plane_2 = float( self.plane_2_entry.get() )
		self.nplanes = float( self.nplanes_entry.get() )
		self.maxhor = float( self.maxhor_entry.get() )
		self.maxhei = float( self.maxhei_entry.get() )
		self.minang = float( self.minang_entry.get() )
		if( self.source_dem in [ 1 , 3 ] ):
			self.vent_lon = float( self.vent_lon_entry.get() )
			self.vent_lat = float( self.vent_lat_entry.get() )
			[ self.pix_matrices , self.pix_matrices_stats , self.pix_matrices_max , self.list_planes ] = pix_height( self.foldername , self.plane_1 , self.plane_2 , self.nplanes , self.maxhor , self.maxhei , self.minang , 1 , self.lon , self.lat , self.vent_lon , self.vent_lat , self.camheight , self.image , self.angles_profile , self.view_profile , self.fix_ang1 , self.fix_ang2 , self.fix_inclination , self.borderx , self.bordery )
		else:
			self.vent_east = float( self.vent_east_entry.get() )
			self.vent_north = float( self.vent_north_entry.get() )
			[ self.pix_matrices , self.pix_matrices_stats , self.pix_matrices_max , self.list_planes ] = pix_height( self.foldername , self.plane_1 , self.plane_2 , self.nplanes , self.maxhor , self.maxhei , self.minang , 2 , self.east , self.north , self.vent_east , self.vent_north , self.camheight , self.image , self.angles_profile , self.view_profile , self.fix_ang1 , self.fix_ang2 , self.fix_inclination , self.borderx , self.bordery )
		self.matrices_available = 1
		self.enabled_disabled()
		messagebox.showinfo( title = None , message = "Matrices created successfully" )

	def plot_pixheight( self ):
		plot_pix_height( self.image , self.pix_matrices , self.pix_matrices_stats , self.pix_matrices_max , self.list_planes )
		
	def enabled_disabled( self ):
		if( self.source_dem_choice == 1 ):
			self.toponame_entry.config( state = 'normal' )
			self.label_toponame.config( state = 'normal' )
		else:
			self.toponame_entry.config( state = 'disabled' )
			self.label_toponame.config( state = 'disabled' )
		if( self.source_dem_choice in [1,3] ):
			self.east_entry.config( state = 'disabled' )
			self.label_east.config( state = 'disabled' )
			self.north_entry.config( state = 'disabled' )
			self.label_north.config( state = 'disabled' )
			self.lon_entry.config( state = 'normal' )
			self.label_lon.config( state = 'normal' )
			self.lat_entry.config( state = 'normal' )
			self.label_lat.config( state = 'normal' )
			self.label_vent_lon.config( state = 'normal' )
			self.vent_lon_entry.config( state = 'normal' )
			self.label_vent_lat.config( state = 'normal' )
			self.vent_lat_entry.config( state = 'normal' )
			self.label_vent_east.config( state = 'disabled' )
			self.vent_east_entry.config( state = 'disabled' )
			self.label_vent_north.config( state = 'disabled' )
			self.vent_north_entry.config( state = 'disabled' )
		else:
			self.lon_entry.config( state = 'disabled' )
			self.label_lon.config( state = 'disabled' )
			self.lat_entry.config( state = 'disabled' )
			self.label_lat.config( state = 'disabled' )
			self.east_entry.config( state = 'normal' )
			self.label_east.config( state = 'normal' )
			self.north_entry.config( state = 'normal' )
			self.label_north.config( state = 'normal' )
			self.label_vent_east.config( state = 'normal' )
			self.vent_east_entry.config( state = 'normal' )
			self.label_vent_north.config( state = 'normal' )
			self.vent_north_entry.config( state = 'normal' )
			self.label_vent_lon.config( state = 'disabled' )
			self.vent_lon_entry.config( state = 'disabled' )
			self.label_vent_lat.config( state = 'disabled' )
			self.vent_lat_entry.config( state = 'disabled' )
		if( self.topo_available == 1 ):
			self.bot_plot_topo.config( state = 'normal' )
			self.bot_create_toptrace.config( state = 'normal' )
		else:
			self.bot_plot_topo.config( state = 'disabled' )
			self.bot_create_toptrace.config( state = 'disabled' )
		if( self.horizon_available == 1 ):
			self.bot_plot_toptrace.config( state = 'normal' )
		else:
			self.bot_plot_toptrace.config( state = 'disabled' )
		if( self.image_available == 1 ):
			self.bot_show_image.config( state = 'normal' )
			self.bot_create_border.config( state = 'normal' )
			self.bot_loadpixheight.config( state = 'normal' )
		else:
			self.bot_show_image.config( state = 'disabled' )
			self.bot_create_border.config( state = 'disabled' )
			self.bot_loadpixheight.config( state = 'disabled' )
		if( self.referenceprofile_available == 1 ):
			self.bot_plot_border.config( state = 'normal' )
		else:
			self.bot_plot_border.config( state = 'disabled' )
		if( self.referenceprofile_available == 1 and self.horizon_available == 1 ):
			self.bot_comp_profile.config( state = 'normal' )
		else:
			self.bot_comp_profile.config( state = 'disabled' )
		if( self.comparison_available == 1 and self.image_available == 1 and self.referenceprofile_available == 1 and self.horizon_available == 1 ):
			self.bot_plot_comp_profile.config( state = 'normal' )
			self.bot_pixheight.config( state = 'normal' )
			self.bot_improve_profile.config( state = 'normal' )
		else:
			self.bot_plot_comp_profile.config( state = 'disabled' )
			self.bot_pixheight.config( state = 'disabled' )
			self.bot_improve_profile.config( state = 'disabled' )
		if( self.matrices_available == 1 ):
			self.bot_plot_pixheight.config( state = 'normal' )
		else:
			self.bot_plot_pixheight.config( state = 'disabled' )

if __name__ == '__main__':
	root = Tk()
	my_gui = MainFrame( root )
	root.after( 100 )
	root.mainloop()
