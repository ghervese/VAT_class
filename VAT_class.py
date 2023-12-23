from astroquery.skyview import SkyView
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy.coordinates import position_angle
import tkinter as tk
from tkinter import ttk
from astropy.wcs import WCS
from astropy.visualization import astropy_mpl_style
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from astropy.visualization.wcsaxes import Quadrangle
from matplotlib.patches import Rectangle
from tkinter import filedialog
import os

#######################################################################################################################################
#                                                      VAT - Virtual Astrophotographer Tool  
#                                    V0.1 by Guillaume Hervé-Secourgeon // herve-guillaume[at]orange.fr
#######################################################################################################################################
# This Python based on astropy and astroquery is dedicated to generate a set of fits images based on NASA's SkyView open observatory
# This class provides a GUI to prepare a set of .fits files that can be post-processed with the prefered softwares of the user.
# In a first frame :
# -----------------
# The user can define the object of interest (ROI) within the Messier, New General Catalog or International Catalog.
# The user can define the field of view (FOV) ot the region of interest in degree.
# On overview is exposed after having selected the object and the radius of the FOV.
# The region of interest is a square area in this first version.
# The metadata are downloaded for the considered ROI.
# In a second frame :
# ------------------
# The user can specify the expected resolution in arcseconds / pixel
# Based on that option a set of tiles is proposed considering a definition of 2400 px X 2400 px for each tile.
# The tiles are superimposed over the overview of the ROI
# The user can also select among the available set of data that are accessible within the diferrent serveys hosted on SkyView server
# The user specify the targeted directory where the data will be uploaded
# The format of the .fits is the following : 
# - For the tile centered on the selected object : 'center_'+survey_name+'_' + object_name + '.fits'
# - For the other tiles : 'tile#_'+survey_name+'_' + object_name + '.fits'
# A short report is exposed, in the prompt, on the screen that summarizes the .fits files that have been downloaded and their location
# on the computer with the the full path.



class VAT_App:
    """
    Classe représentant une application d'aperçu des objets célestes.

    Cette application utilise tkinter pour créer une interface graphique permettant à l'utilisateur de spécifier
    les paramètres d'un objet céleste, de récupérer les données correspondantes et d'afficher un aperçu de l'objet.
    """
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("VAT - Virtual Astrophographer Tool by G. Hervé-Secourgeon")
        self.fig_photo = None
        self.fig_photo_tiles = None
        self.nb_pixels = 2000

    def retrieve_data(self):
        """
        Récupère les données saisies par l'utilisateur dans les zones de texte.

        Cette méthode récupère les valeurs de la longueur focale, de la taille du capteur, de la taille du pixel
        et du nom de l'objet à partir des zones de texte de l'interface graphique. Si aucune valeur n'est saisie,
        des valeurs par défaut sont utilisées. Les données récupérées sont ensuite imprimées dans la console.
        """
        self.angle_obs = float(self.angle_obs_entry.get()) if self.angle_obs_entry.get() else 2.
        self.fov = self.angle_obs*u.deg
        self.object = self.object_entry.get() if self.object_entry.get() else "M31"
        catalogues_supportes = ['NGC', 'M', 'IC']
        catalogue = None
        for c in catalogues_supportes:
            if self.object.startswith(c):
                catalogue = c
                break

        if catalogue is None:
            print("Catalogue non supporté. Veuillez utiliser NGC, M ou IC.")
            return
    
    def retrieve_survey(self):
        """
        Récupère les données saisies par l'utilisateur dans les zones de texte.

        Cette méthode récupère les valeurs de la longueur focale, de la taille du capteur, de la taille du pixel
        et du nom de l'objet à partir des zones de texte de l'interface graphique. Si aucune valeur n'est saisie,
        des valeurs par défaut sont utilisées. Les données récupérées sont ensuite imprimées dans la console.
        """
        self.survey = self.object_entry.get() if self.object_entry.get() else "DSS2 Blue"
        available_surveys = ['DDS1 Blue', 'DSS2 Blue', 'DSS2 Red', 'DSS2 IR']
        survey_catalogue = None
        for c in available_surveys:
            if self.survey.startswith(c):
                survey_catalogue = c
                break

        if survey_catalogue is None:
            print("That survey is not available. Please use : DDS1 Blue, DSS2 Blue, DSS2 Red', DSS2 IR.")
            return

    def SNR_calculation(self, picture):
        """
        Calcule le rapport signal/bruit d'une image et le niveau de bruit moyen.

        Cette méthode prend en entrée une image sous forme de tableau numpy et calcule le rapport signal/bruit en
        utilisant la moyenne du signal et l'écart type du bruit. Elle retourne également le niveau de bruit moyen.

        Args:
            picture (ndarray): Tableau numpy représentant une image.

        Returns:
            tuple: Rapport signal/bruit et niveau de bruit moyen.
        """
        signal_mean = np.mean(picture)
        noise = picture - signal_mean
        noise_mean = np.mean(noise)
        noise_std = np.std(noise)
        snr = signal_mean / noise_std
        return snr, noise_mean

    def calculate_number_tiles(self):
        """
        Calcul le nombre de tuiles en considérant la taille de la ROI, la résolution cible et le nombre de pixels par tuile

        Args:
            nb_pixels : nombre de pixels par tuile
            self.angle_obs : valeur de l'angle de la ROI en °
            self.tile_resol : valeur de la résolution dans une tuile en secondes d'arc par pixels

        Returns:
            self.tile_fov : l'angle de vue d'une tuile élémentaire
            self.number_tiles : le nombre de tuiles par ROI dans chaque direction => le nombre total de tuiles est au carré
        """
        nb_pixels = self.nb_pixels
        self.tile_resol = float(self.resolution_entry.get()) if self.resolution_entry.get() else 0.7
        
        self.cover_pct = float(self.cover_entry.get()) if self.resolution_entry.get() else 10.
        # Conversion de la résolution de "/px en °/px
        self.tile_resol_deg = self.tile_resol / 3600.
        # Calcul de la taille (FOV_tile) d'une tuile élémentaire pour satisfaire le critère de résolution
        self.tile_fov = self.tile_resol / 3600. * nb_pixels
        print("The FOV of global target has been set to :"+str(self.angle_obs)+" in °")
        print("The targeted resolution per tile has been set by user to :"+str(self.tile_resol)+" in \"/px")
        print("The FOV of the tile in order to satisfy the resolution criterion is set to: "+str(np.round(self.tile_fov,2))+" in °")
        # Calcul du nombre de tuiles de la FOV élémentaire déterminée pour couvrir la FOV de la ROI de l'objet choisi
        self.number_tiles = np.int64(np.ceil(self.angle_obs / self.tile_fov - self.cover_pct/100.) * 1./(1. - self.cover_pct/100.))
        print("The number of tiles necessary to comply with the FOV of the ROI with a cover of "+str((self.cover_pct))+"% and 2400x2400 px² is set to :"+str(self.number_tiles**2))
        print("The number of tiles per direction is then set to "+str(self.number_tiles))

    def tile_coordinates(self):
        """
        Calcul en utilisant, la méthode directional_offset_by, des coordonnées des différentes tuiles pour :
            - prévisualiser leur emplacement sur l'objet à imager
            - préparer l'obtention des .fits avec SkyView

        Args:
            self.object: Nom de l'objet ciblé
            self.number_tiles : Nombre de tuiles
            self.tile_fov : Angle de vue de chaque tuile
        Returns:

        """
        print("Execution de tile coordinates")
        self.center_coords = SkyCoord.from_name(self.object)
        self.offset_value = (self.tile_fov*(1-self.cover_pct/100.) )* u.deg
        # Direction de référence pour le décalage suivant les 
        # déclinaisons positives, on a pris pi/4 arbitrairement. Il fallait logiquement une valeur entre 0. et pi/2.
        # Le même raisonnement est appliquée pour chaque direction de décalage. Seul le signe change et les coordonnées du 2nd
        # point, on inverse simplement longitude et latitude.
        position_angle_DECplus = position_angle(0.,0.,0.,np.pi/4.) 
        position_angle_DECmoins = position_angle(0.,0.,0.,-np.pi/4.)
        position_angle_RAplus = position_angle(0.,0.,np.pi/4.,0.)
        position_angle_RAmoins = position_angle(0.,0.,-np.pi/4.,0.)
        # Translation du point de départ à partir des corrdonnées du centre de l'objet ciblé
        coord_starting_point1 = self.center_coords.directional_offset_by(position_angle_RAmoins,self.offset_value/2.*float(self.number_tiles))
        coord_starting_point = coord_starting_point1.directional_offset_by(position_angle_DECmoins,self.offset_value/2.*float(self.number_tiles))
        # On applique ensuite le décalage pour déterminer l'emplacement de chaque tuile à partir de cette nouvelle origine
        self.tile_coordinates_center = []
        for i in range(self.number_tiles):
            for j in range(self.number_tiles):
                coord_tile_step1 = coord_starting_point.directional_offset_by(position_angle_RAplus,self.offset_value*float(i))
                coord_tile_final = coord_tile_step1.directional_offset_by(position_angle_DECplus,self.offset_value*float(j))
                self.tile_coordinates_center.append(coord_tile_final)
                print("Coordinates of tile "+str(i+1)+"x"+str(j+1)+" :"+str(coord_tile_final))
        print(self.tile_coordinates_center)
        print(np.shape(self.tile_coordinates_center))

    def data_overview(self):
        """
        Génère un aperçu de l'objet céleste en utilisant les données récupérées.

        Cette méthode vérifie le catalogue (NGC, M, IC) auquel appartient l'objet céleste spécifié. Si le catalogue
        est valide, elle récupère les informations sur l'objet à partir du catalogue SkyView, et récupère les images
        correspondantes à partir du service SkyView. Elle effectue également des calculs de rapport signal/bruit.
        Enfin, elle affiche l'aperçu de l'objet dans l'interface graphique.
        """
        nb_pixels = self.nb_pixels # This parameter is only adjustable through the edition of class. It is not an option accessible with the GUI

        if self.object is not None:
            # images = SkyView.get_images(self.object, survey='DSS')[0][0]
            # image_Ha = SkyView.get_images(self.objet, survey='H-Alpha Comp')[0][0]
            hdu = SkyView.get_images(self.object, survey=['DSS'],pixels=np.int64(np.round(nb_pixels/2.)),radius=self.fov)[0][0]
            wcs = WCS(hdu.header)
            snr_visible, noise_mean_visible = self.SNR_calculation(hdu.data)

            # print("Rapport signal/bruit H-alpha ", self.objet, " :", snr_Ha)
            # print("Niveau de bruit H-alpha de ", self.objet, " :", noise_mean_Ha)
            print("Noise to signal ratio of the overview ", self.object, " :", snr_visible)
            print("Mean noise value of the overview   ", self.object, " :", noise_mean_visible)

            plt.clf()
            fig = plt.figure()
            plt.style.use(astropy_mpl_style)
#            ax = fig.add_subplot(projection=wcs, label='overlays')
            ax = fig.add_subplot(projection=wcs)
            
            ax.imshow(hdu.data, cmap='binary_r', origin='lower')


            """
            Lis of available color maps:
            'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 
            'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 
            'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 
            'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 
            'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 
            'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 
            'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 
            'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 
            'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 
            'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 
            'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 
            'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 
            'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 
            'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 
            'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 
            'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 
            'viridis_r', 'winter', 'winter_r'
            """

            if self.fig_photo is not None:
                self.fig_photo.get_tk_widget().destroy()
            self.fig_photo = FigureCanvasTkAgg(fig, master=self.canvas)
            self.fig_photo.draw()
            self.fig_photo.get_tk_widget().pack()

        #    hdu.writeto('DSS_' + self.object + '.fits', overwrite=True)
        else:
            print("Objet non trouvé dans le catalogue.")


    def data_overview_tiles(self):
        """
        Génère un aperçu de l'objet céleste en utilisant les données récupérées.

        Cette méthode vérifie le catalogue (NGC, M, IC) auquel appartient l'objet céleste spécifié. Si le catalogue
        est valide, elle récupère les informations sur l'objet à partir du catalogue Simbad, calcule les paramètres
        nécessaires pour l'aperçu (résolution, taille de l'image, etc.) et récupère les images correspondantes à
        partir du service SkyView. Elle effectue également des calculs de rapport signal/bruit. Enfin, elle affiche
        l'aperçu de l'objet dans l'interface graphique.
        """
        nb_pixels = self.nb_pixels # This parameter is only adjustable through the edition of class. It is not an option accessible with the GUI

        if self.object is not None:
            # images = SkyView.get_images(self.object, survey='DSS')[0][0]
            # image_Ha = SkyView.get_images(self.objet, survey='H-Alpha Comp')[0][0]
            hdu = SkyView.get_images(self.object, survey=['DSS'],pixels=np.int64(np.round(nb_pixels/2.)),radius=self.fov)[0][0]
            wcs = WCS(hdu.header)
            print(wcs)
            snr_visible, noise_mean_visible = self.SNR_calculation(hdu.data)

            # print("Rapport signal/bruit H-alpha ", self.objet, " :", snr_Ha)
            # print("Niveau de bruit H-alpha de ", self.objet, " :", noise_mean_Ha)
            print("Noise to signal ratio of the overview ", self.object, " :", snr_visible)
            print("Mean noise value of the overview   ", self.object, " :", noise_mean_visible)

            plt.clf()
            fig_tiles = plt.figure()
            plt.style.use(astropy_mpl_style)
            #ax = fig_tiles.add_subplot(projection=wcs, label='overlays')
            ax = fig_tiles.add_subplot(projection=wcs)
            print(self.tile_fov)
            ax.imshow(hdu.data, cmap='binary_r', origin='lower')
            ax.grid(color='black',ls='solid')
            for i in range(len(self.tile_coordinates_center)):                
                r = Quadrangle(tuple([self.tile_coordinates_center[i].ra.degree,self.tile_coordinates_center[i].dec.degree])*u.deg, 1.4*self.tile_fov*u.deg, self.tile_fov*u.deg,
#                r = Rectangle(tuple([self.tile_coordinates_center[i].ra.degree,self.tile_coordinates_center[i].dec.degree]), self.tile_fov*1.4, self.tile_fov,
                edgecolor='red', facecolor='none',
                transform=ax.get_transform('world')
                )
                ax.add_patch(r)



            if self.fig_photo_tiles is not None:
                self.fig_photo_tiles.get_tk_widget().destroy()
            self.fig_photo_tiles = FigureCanvasTkAgg(fig_tiles, master=self.canvas_tiles)
            self.fig_photo_tiles.draw()
            self.fig_photo_tiles.get_tk_widget().pack()

        #    hdu.writeto('DSS_' + self.object + '.fits', overwrite=True)
        else:
            print("Objet non trouvé dans le catalogue.")


    def select_menu_fit1(self,event):
        self.survey_fit1 = self.list_combo_survey_1.get()
    def select_menu_fit2(self,event):
        self.survey_fit2 = self.list_combo_survey_2.get()
    def select_menu_fit3(self,event):
        self.survey_fit3 = self.list_combo_survey_3.get()
    def select_menu_fit4(self,event):
        self.survey_fit4 = self.list_combo_survey_4.get()

    def create_subfolder(self):
        source_path = filedialog.askdirectory(title='Select the Parent Directory')
        path = os.path.join(source_path, 'Images')
        os.makedirs(path)

    def import_data_from_SkyView(self):
        for i in range(len(self.tile_coordinates_center)):
            nb_pixels = self.nb_pixels
            image_channel_1 = SkyView.get_images(self.tile_coordinates_center[i], survey=self.survey_fit1,pixels=nb_pixels,radius=self.tile_fov*u.deg)[0][0]
            image_channel_2 = SkyView.get_images(self.tile_coordinates_center[i], survey=self.survey_fit2,pixels=nb_pixels,radius=self.tile_fov*u.deg)[0][0]
            image_channel_3 = SkyView.get_images(self.tile_coordinates_center[i], survey=self.survey_fit3,pixels=nb_pixels,radius=self.tile_fov*u.deg)[0][0]
            image_channel_4 = SkyView.get_images(self.tile_coordinates_center[i], survey=self.survey_fit4,pixels=nb_pixels,radius=self.tile_fov*u.deg)[0][0]
            image_channel_1.writeto(self.object + '_' + self.survey_fit1 + '_tile_' + str(i)+ '.fits', overwrite=True)
            image_channel_2.writeto(self.object + '_' + self.survey_fit2 + '_tile_' + str(i)+ '.fits', overwrite=True)
            image_channel_3.writeto(self.object + '_' + self.survey_fit3 + '_tile_' + str(i)+ '.fits', overwrite=True)
            image_channel_4.writeto(self.object + '_' + self.survey_fit4 + '_tile_' + str(i)+ '.fits', overwrite=True)

    def create_gui(self):
        """
        Crée l'interface graphique de l'application.

        Cette méthode utilise la bibliothèque tkinter pour créer les éléments graphiques de l'interface utilisateur,
        tels que les zones de texte, les boutons et les onglets. Elle associe également les méthodes appropriées
        aux événements des boutons.
        """
        # Création de l'onglet principal
        tab_control = ttk.Notebook(self.root)
        tab1 = ttk.Frame(tab_control)
        tab_control.add(tab1, text='Target overview')
        tab_control.pack(expand=1, fill='both')

        # Création de l'onglet de viusalisation du rapport signal sur bruit des différents filtres
        tab2 = ttk.Frame(tab_control)
        tab_control.add(tab2, text='FITS content specification')
        tab_control.pack(expand=1, fill='both')

        # Création de l'onglet de récupération des données
        tab3 = ttk.Frame(tab_control)
        tab_control.add(tab3, text='Download data')
        tab_control.pack(expand=1, fill='both')



        ##########################################################
        # Code relatif au premier onglet
        ##########################################################

        # Cadre pour les entrées utilisateur
        input_frame = ttk.LabelFrame(tab1, text='Region of Interest')
        input_frame.pack(expand=1, fill='both', padx=2, pady=2)

        # Zone de texte pour le nom de l'objet
        object_label = ttk.Label(input_frame, text='Targeted object (NGC/M/IC #):')
        object_label.grid(row=0, column=0, padx=5, pady=5)
        self.object_entry = ttk.Entry(input_frame)
        self.object_entry.insert(0, 'M31')
        self.object_entry.grid(row=0, column=1, padx=5, pady=5)

        # Zone de texte pour la longueur focale
        angle_obs_label = ttk.Label(input_frame, text='Fielf of vision in °:')
        angle_obs_label.grid(row=1, column=0, padx=5, pady=5)
        self.angle_obs_entry = ttk.Entry(input_frame)
        self.angle_obs_entry.insert(0, '2.')
        self.angle_obs_entry.grid(row=1, column=1, padx=5, pady=5)

        # Bouton pour récupérer les données
        retrieve_button = ttk.Button(input_frame, text='Retrieve data', command=self.retrieve_data)
        retrieve_button.grid(row=0, column=2, columnspan=2, padx=5, pady=5)

        # Bouton pour générer l'aperçu
        preview_button = ttk.Button(input_frame, text='Generate overview', command=self.data_overview)
        preview_button.grid(row=1,column=2,columnspan=2, padx=5, pady=5)

        # Bouton pour quitter
        exit_button = ttk.Button(input_frame, text="Quit", command=self.Close)
        exit_button.grid(row=0,column=4,columnspan=2,padx=5,pady=5)

        # Cadre pour l'aperçu de l'objet
        preview_frame = ttk.LabelFrame(tab1, text='Overview of the global Region of Interest')
        preview_frame.pack(expand=1, fill='both', padx=2, pady=2)

        # Canvas pour afficher l'aperçu
        self.canvas = tk.Canvas(preview_frame)
        self.canvas.pack(expand=1, fill='both')

        ######################################################
        # Code pour le deuxième onglet
        ######################################################
        # Cadre pour les paramètres de tuiles
        input_tiles = ttk.LabelFrame(tab2, text='Parameters of the expected tiles')
        input_tiles.pack(expand=1, fill='both', padx=2, pady=2)
        # Cadre pour visualiser les résultats dans un graphique
        view_tiles = ttk.LabelFrame(tab2, text='Overview of the locations of the tiles on the object')
        view_tiles.pack(expand=1, fill='both', padx=2, pady=2)
        # Zone de texte pour la resolution qui pilotera la taille des tuiles en considérant une taille de 2400x2400 px²
        resolution_label = ttk.Label(input_tiles, text='Resolution in "/px :')
        resolution_label.grid(row=0, column=0, padx=5, pady=5)
        self.resolution_entry = ttk.Entry(input_tiles)
        self.resolution_entry.insert(0, '0.7')
        self.resolution_entry.grid(row=0, column=1, padx=5, pady=5)
        # Pourcentage de recouvrement entre 2 tuiles
        cover_label = ttk.Label(input_tiles, text='Cover percentage between tiles:')
        cover_label.grid(row=0, column=2, padx=5, pady=5)
        self.cover_entry = ttk.Entry(input_tiles)
        self.cover_entry.insert(0, '10.')
        self.cover_entry.grid(row=0, column=3, padx=5, pady=5)
        # Bouton pour calculer le nombre de tuiles
        tiles_number_button = ttk.Button(input_tiles, text='Calculate tiles', command=self.calculate_number_tiles)
        tiles_number_button.grid(row=1, column=0, columnspan=2, padx=5, pady=5)

        # Bouton pour calculer les coordonnées des tuiles sur la ROI de la cible
        tiles_preview_button = ttk.Button(input_tiles, text='Preview of the tiles', command=self.tile_coordinates)
        tiles_preview_button.grid(row=1, column=1, columnspan=2, padx=5, pady=5)

        # Bouton pour visualiser les tuiles sur la ROI de la cible
        tiles_preview_button = ttk.Button(input_tiles, text='Plot tiles over preview', command=self.data_overview_tiles)
        tiles_preview_button.grid(row=1, column=2, columnspan=2, padx=5, pady=5)


        



        # Canvas pour afficher l'aperçu
        self.canvas_tiles = tk.Canvas(view_tiles)
        self.canvas_tiles.pack(expand=1, fill='both')


        ######################################################
        # Code pour le troisième onglet
        ######################################################
        # Cadre pour les longueurs d'onde / filtres
        input_surveys = ttk.LabelFrame(tab3, text='Parameters of the surveys')
        input_surveys.pack(expand=1, fill='both', padx=2, pady=2)
        





        # Zone de texte pour la taille du pixel
        dep1_label = ttk.Label(input_surveys, text='Survey for Channel 1:')
        dep1_label.grid(row=0, column=0, padx=5, pady=5)
        # liste des unites
        self.list_survey_fit = ['DSS','DSS1 Blue','DSS2 Blue','DSS1 Red','DSS2 Red','DSS2 IR']
        # creation comboBox
        self.list_combo_survey_1=ttk.Combobox(input_surveys, values=self.list_survey_fit)
        self.list_combo_survey_1.current(0)
        #Position de la ComboBox
        self.list_combo_survey_1.grid(row=1, column=0, padx=5, pady=5)
        self.list_combo_survey_1.bind("<<ComboboxSelected>>",self.select_menu_fit1)


        # Zone de texte pour la taille du pixel
        dep1_label = ttk.Label(input_surveys, text='Survey for Channel 2:')
        dep1_label.grid(row=0, column=1, padx=5, pady=5)
        # liste des unites
        self.list_survey_fit = ['DSS','DSS1 Blue','DSS2 Blue','DSS1 Red','DSS2 Red','DSS2 IR']
        # creation comboBox
        self.list_combo_survey_2=ttk.Combobox(input_surveys, values=self.list_survey_fit)
        self.list_combo_survey_2.current(2)
        #Position de la ComboBox
        self.list_combo_survey_2.grid(row=1, column=1, padx=5, pady=5)
        self.list_combo_survey_2.bind("<<ComboboxSelected>>",self.select_menu_fit2)

       # Zone de texte pour la taille du pixel
        dep1_label = ttk.Label(input_surveys, text='Survey for Channel 3:')
        dep1_label.grid(row=0, column=2, padx=5, pady=5)
        # liste des unites
        self.list_survey_fit = ['DSS','DSS1 Blue','DSS2 Blue','DSS1 Red','DSS2 Red','DSS2 IR']
        # creation comboBox
        self.list_combo_survey_3=ttk.Combobox(input_surveys, values=self.list_survey_fit)
        self.list_combo_survey_3.current(4)
        #Position de la ComboBox
        self.list_combo_survey_3.grid(row=1, column=2, padx=5, pady=5)
        self.list_combo_survey_3.bind("<<ComboboxSelected>>",self.select_menu_fit3)

       # Zone de texte pour la taille du pixel
        dep1_label = ttk.Label(input_surveys, text='Survey for Channel 4:')
        dep1_label.grid(row=0, column=3, padx=5, pady=5)
        # liste des unites
        self.list_survey_fit = ['DSS','DSS1 Blue','DSS2 Blue','DSS1 Red','DSS2 Red','DSS2 IR']
        # creation comboBox
        self.list_combo_survey_4=ttk.Combobox(input_surveys, values=self.list_survey_fit)
        self.list_combo_survey_4.current(5)
        #Position de la ComboBox
        self.list_combo_survey_4.grid(row=1, column=3, padx=5, pady=5)
        self.list_combo_survey_4.bind("<<ComboboxSelected>>",self.select_menu_fit4)

        # Attribution des valeurs pour défaut au cas où l'utilisateur souhaite garder les paramètres proposés par défaut
        self.survey_fit1 = self.list_survey_fit[0]
        self.survey_fit2 = self.list_survey_fit[2]
        self.survey_fit3 = self.list_survey_fit[4]
        self.survey_fit4 = self.list_survey_fit[5]



        chose_path_button = ttk.Button(input_surveys, text="Select a Folder", command=self.create_subfolder)
        chose_path_button.grid(row=2, column=0, columnspan=2, padx=5, pady=5)

        launch_import_button = ttk.Button(input_surveys, text="Import .fits Data", command=self.import_data_from_SkyView)
        launch_import_button.grid(row=2, column=2, columnspan=2, padx=5, pady=5)

    def Close(self):
        self.root.destroy()
    
    def run(self):
        """
    Lance l'application.

        Cette méthode démarre la boucle principale de l'interface graphique tkinter et affiche l'application
        à l'écran.
        """
        self.create_gui()
        self.root.mainloop()
# Instanciation de l'application et exécution
app = VAT_App()
app.run()