import csv
import os
import sys
import numpy as np
from mirpy import miriad
import pyfits
import matplotlib
import matplotlib.pyplot as plt


# ---------------------------------

def get_tsys(kw, current_time, ska_par):
    # a function to get the system temperature:
    if kw['tel'] == 'SKA1_Mid':
        eff = float(ska_par['mid_eff'])
        dish = float(ska_par['mid_dish'])
        if kw['band'] == 'band1':
            tsys_file = str(ska_par['mid_tsys_band1'])
        if kw['band'] == 'band2':
            tsys_file = str(ska_par['mid_tsys_band2'])
           
    if kw['tel'] == 'SKA1_Low':
        eff = float(ska_par['low_eff'])
        dish = float(ska_par['low_dish'])
        tsys_file = str(ska_par['low_tsys_band1'])
    
    #os.system('less calculator/data/mid_tsys_band2.csv > less_temp.txt')
    f=open(tsys_file,"rb")
    #f=open("'calculator/data/mid_tsys_band2.csv'",'rb')
    reader = csv.reader(f)
    tsys = list(reader)
    print tsys

    if kw['mode'] == 'line':
        freq = float(kw['freq'])
        for i in range(1,len(tsys)-1,1):
            
            if freq >= float(tsys[i][0]) and freq <= float(tsys[i+1][0]):
                if freq == float(tsys[i][0]):
                    A_tsys = float(tsys[i][1])
                else:
                    A_tsys = (float(tsys[i][1]) + float(tsys[i+1][1])) / 2 
                                                        
    if kw['mode'] == 'cont':
        freq1 = float(kw['freq']) - float(kw['width'])/2
        freq2 =float(kw['freq']) + float(kw['width'])/2
        for i in range(1,len(tsys)-1,1):
            if freq1 >= float(tsys[i][0]) and freq1 <= float(tsys[i+1][0]):
                i1 = i
            if freq2 >=float(tsys[i][0]) and freq2 <= float(tsys[i+1][0]):
                i2 = i
                
        A_tsys = 0
        if i1 == i2:
            A_tsys = tsys[i1][1]
        else:
            for i in range(i1,i2,1):
                A_tsys = A_tsys + float(tsys[i][1])
            A_tsys = A_tsys / ( i2-i1 + 1.)

    output = np.pi*(dish/2)**2 * eff / A_tsys
    #output = np.pi*(dish/2)**2 / A_tsys                              
    return output

#___________________________________________________    

def make_uv(kw, current_time, ska_par):
    # A function to generate simulated visibility data
    #os.system('touch bla_' + current_time + '.txt')
    # make a directory for this run:
    outdir = '/tmp/ska_calculator/data/runs/' + current_time + '/'
    #os.system('mkdir ' + outdir)
    #outdir = '/tmp/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outname = outdir + 'data.uv'

    # select the telescope configuration
    if kw['tel'] == 'SKA1_Mid':
        config_file = ska_par['mid_config']
        baseunit = str(ska_par['mid_baseunit'])
        dish = float(ska_par['mid_dish'])
        eff = str(ska_par['mid_eff'])
        #tsys = str(ska_par['mid_tsys'])
        lat = str(ska_par['mid_lat'])
        eff = float(ska_par['mid_eff'])
        sys_eff = float(ska_par['mid_sys_eff'])

    if kw['tel'] == 'SKA1_Low':
        config_file = ska_par['low_config']
        baseunit = str(ska_par['low_baseunit'])
        dish = float(ska_par['low_dish'])
        eff = str(ska_par['low_eff'])
        #tsys = str(ska_par['mid_tsys'])
        lat = str(ska_par['low_lat'])
        eff = float(ska_par['low_eff'])
        sys_eff = float(ska_par['low_sys_eff'])

    

    tsys = get_tsys(kw, current_time, ska_par)     

    #AP!! This step can be faster by generating a pre-defined config file    
    #read the .csv file with the telescope configuration    
    #f  = open(config_file, "rb")
    #reader = csv.reader(f)
    #configuration = list(reader)
    #f.close()    

    #f = open('calculator/data/configurations/blaa.cfg','w')
    #for i in range(len(configuration)):
    #    # change the order of x and y, as this is the default for miriad
    #    #f.write(configuration[i][0] + '   ' + configuration[i][1]  + '\n')
    #    f.write(configuration[i][1] + '   ' + configuration[i][0]  + '\n')
    #f.flush()
    #os.fsync(f.fileno())
    #f.close()
    
    # calculate the cycle time based on the taper:
    # Theta = lambda / D
    # D = lambda / Theta
    # D = (0.3/freq) / (Theta''/3600 * pi/180)
    # D = (0.3 * 3600 *180) / (Theta'' * freq * pi)
    # D = 61879/(Theta'' * freq)
    # D approximates the longest baseline, relevant for the chosen baseline
    # calculate the time, in which the array moves 15 meters on this longes baseline
    # v_tel = (2*pi*D/2) / (24*3600) [m/s]
    # to move 7.5m (15m dish and Nyquist sampling):
    # t = 7.5/v_tel = 206265/D
    # t = (Theta'' * freq) * 3.33
        
    # !! Double check this and look at actual differences in results    

    #cycle_time = 60# seconds
    #cycle_time = float(kw['taper']) * float(kw['freq']) * 3.33    
    #cycle_time =704
    
    cycle_time = ska_par['sample_time']    


    cycle_time = float(cycle_time) / 3600 #hours
    frequency = float(kw['freq'])/1000 # convert MHz to GHz
    aperture = np.pi * (dish/2.0)**2
    jyperk_value = (2*1.38e3)/(eff * aperture * sys_eff)
    empty_file='ska_calculator/data/empty.dat'

    # if the name of the configuration file is too long, miriad crashes
    os.system('cp ' + config_file + ' /tmp/temp.cfg')



    if kw['mode'] == 'line':
        corr_setup = '1,1,0,' + kw['width']
    if kw['mode'] == 'cont':
        corr_setup = str(ska_par['cont_chan']) + ',1,0,' + kw['width']

    miriad.uvgen(source=empty_file,
        ant='/tmp/temp.cfg',
        baseunit=str(baseunit),
        #telescop='equatorial,0',
        telescope='altaz,0',
        lat=str(lat),
        corr=corr_setup,
        freq=str(frequency),
        radec='0,' + str(kw['dec']),
        harange=str(kw['hamin']) + ',' + str(kw['hamax']) + ',' + str(cycle_time),
        systemp=str(tsys),
        jyperk=str(jyperk_value),
        out=outname)

    #miriad.uvdump(vis=outname,
    #    vars='uu,vv',
    #    log='/tmp/bla.log')

    #os.system('invert vis=' + outname + ' map=bla.map imsize=512 cellsize=1- options=mfs robust=-2')

    #miriad.invert(vis=outname,
    #    map = '/tmp/bla.map',
    #    imsize=512,
    #    cellsize=10,
    #    options='mfs',
    #    robust=-2)


    #os.system('rm -rf /tmp/temp.cfg')
    return 

#______________________________________________

def make_image(kw,current_time,ska_par):
    # function to make a noise map and psf
    outdir = '/tmp/ska_calculator/data/runs/' + current_time + '/'

    uvfile = outdir + 'data.uv'
    map = outdir + 'noise.im'
    psf = outdir + 'psf.im'

    os.system('chmod -R a+w ' + outdir)
    os.system('chmod a+x ' + uvfile)

    
    print 'uvfile: ' + uvfile
    print 'map: ' + map
    print 'psf: ' + psf


    if kw['tel'] == 'SKA1_Mid':
        dish = float(ska_par['mid_dish'])
        
    if kw['tel'] == 'SKA1_Low':
        dish = float(ska_par['low_dish'])
        # do something else
        

    #pixels=2048
    # calculate the number of pixels based on the resolution and cell size
    taper = float(kw['taper'])
    cellsize = taper/4 
    
    # FoV = lambda/D * 1/1.133, where lambda = 0.3/freq and D is the telescope diameter
    FoV = 300 / (1.133 * float(kw['freq']) * dish) * 180 / np.pi
    npix = FoV / (cellsize / 3600)
    
    # calculate the cycle time for perfect sampling
    # D = 182/taper ; diameter in  kilo lamba
    D = 182./taper * 1000 * 300/float(kw['freq']) # Dimater in meters
    tel_speed = np.pi * D / (24*3600) # telescope speed in m/s
    cycle_time = dish/tel_speed
    if cycle_time < float(ska_par['sample_time']):
        time_ratio = float(ska_par['sample_time']) / cycle_time
    else:
        time_ratio = 1

    scaling_factor = time_ratio**0.25
    pix_factor = 2**scaling_factor
    if npix > 64:      
        scale_pix = np.ceil(npix /  pix_factor)
    else:
        scale_pix = npix
                                     

    robust=-2 # uniform
    #robust=2 # natural
    #robust=0
    

    # lambda/D = Theta
    # lambda/D = taper/3600 * pi/180
    # D =  3600 * 180 lambda / taper /pi
    # lambda = 0.3/freq
    # D = 61879/(taper * freq)

    # klambda = D/(1000*lambda) = (D * freq) /(300)
    # klambda = 206.26/taper

    # klambda = 182/taper # see miriad documentation of fwhm in invert

    klambda = 182 / taper 
    uvrange = 'uvrange(0,' + str(klambda) + ')'

    if kw['mode'] == 'line':
        line_setup = 'channel,1,1,1,1'
        options_setup = ''
    if kw['mode'] == 'cont':
        line_setup = 'channel,1,1,' + str(ska_par['cont_chan']) + ',1'
        options_setup = 'mfs'

    #miriad.invert(vis=uvfile,
    #             map=map,
    #             beam=psf,
    #             #imsize=pixels,
    #             imsize=scale_pix,
    #             cell=cellsize,
    #             robust=robust,
    #             fwhm=taper,
    #             #select=uvrange,
    #             line=line_setup,
    #             options=options_setup)

    invert_command = 'invert vis=' + uvfile + ' map=' + map + ' beam=' + psf + ' imsize=' + str(scale_pix) + ' cell=' +str(cellsize) + ' robust=' + str(robust) + ' fwhm=' + str(taper) + ' line=' + line_setup + ' options=' + options_setup
    print invert_command


    os.system(invert_command)
            
    return


def plot_psf(kw,current_time):
    # store a plot of the psf
    outdir = '/tmp/ska_calculator/data/runs/' + current_time + '/'
    psf = outdir + 'psf.im'
    #public_dir = 'ska_calculator/public/data/' + current_time + '/'
    public_dir = '/home/ec2-user/turbogears/ska-calculator/ska_calculator/public/data/' + current_time + '/'
    psf_fits = outdir + 'psf.fits'
    plotname = outdir + '/psf.png'
    
    
    if not os.path.exists(public_dir):
            os.makedirs(public_dir)


    print 'psf=', psf
    print 'public_dir=', public_dir

    miriad.fits(In=psf, out=psf_fits, op='xyout')
    f=pyfits.open(psf_fits)
    beam_data = f[0].data
    beam_data = beam_data[0,:,:]
    #beam_stat = beam_data[960:-960,960:-960]
    beam_stat = beam_data
    min_range = beam_stat.min()
    cdelt = f[0].header['cdelt1']*3600
    f.close()
    
    plt.figure(1, figsize=(8,6))
    plt.clf()  # clear the previous figure
    
    plt.subplot(1,1,1)
    midpoint = np.round(len(beam_data)/2)
    box_rad = 16
    box_min = midpoint-box_rad
    box_max = midpoint+box_rad
    extent_min = -box_rad * cdelt
    extent_max = box_rad * cdelt
    #extent = [-64*cdelt,64*cdelt,-64*cdelt,64*cdelt]
    #plt.imshow(beam_data[960:-960,960:-960], extent=extent,interpolation='nearest',vmin=min_range,vmax=-3*min_range)
    extent = [extent_min,extent_max,extent_min,extent_max]
    plt.imshow(beam_data[box_min:box_max,box_min:box_max], extent=extent,interpolation='nearest',vmin=min_range,vmax=-3*min_range)
    plt.colorbar()

    levels1=[0.1,0.2,0.5]
    levels2=[-0.05, -0.01]
    plt.contour(beam_data[box_min:box_max,box_min:box_max], levels1, hold='on', colors = 'k',origin='upper', extent=extent, aspect='auto')
    plt.contour(beam_data[box_min:box_max,box_min:box_max], levels2, hold='on', colors = 'w',origin='upper', extent=extent, aspect='auto')
    plt.xlabel('arcsec', size=10)
    plt.ylabel('arcsec', size=10)
    plt.title('Beam PSF', size=12)
    plt.savefig(plotname)

    #print 'BLAAAAAAAAAAAAAAAA$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'

    #print plotname
    #print public_dir
    #os.system('ls ' + public_dir)
    
    os.system('cp ' + plotname + ' ' + public_dir + '/.')
    os.system('cp ' + psf_fits + ' ' + public_dir + '/.')
    
    #os.system('ls -l ' + public_dir)

    return

#_______________________________________________

def plot_noise(kw,current_time):
    # store a plot of the psf
    outdir = '/tmp/ska_calculator/data/runs/' + current_time + '/'
    data = outdir + 'noise.im'
    data_fits = outdir + 'noise.fits'
    plotname = outdir + 'noise.png'

    #public_dir = 'ska_calculator/public/data/' + current_time + '/'
    public_dir = '/home/ec2-user/turbogears/ska-calculator/ska_calculator/public/data/' + current_time +'/'

    miriad.fits(In=data, out=data_fits, op='xyout')
    f=pyfits.open(data_fits)
    noise_data = f[0].data
    noise_data = noise_data[0,0,:,:]
    cdelt = f[0].header['cdelt1']*3600
    f.close()

    plt.figure(1, figsize=(8,6))
    plt.clf()  # clear the previous figure
    plt.subplot(1,1,1)
    midpoint = np.round(len(noise_data)/2)
    box_rad = 64
    box_min = midpoint-box_rad
    box_max = midpoint+box_rad
    extent_min = -box_rad * cdelt
    extent_max = box_rad * cdelt
    extent = [extent_min,extent_max,extent_min,extent_max]
    plt.imshow(noise_data[box_min:box_max,box_min:box_max], extent=extent,interpolation='nearest')
    plt.colorbar()
    plt.xlabel('arcsec', size=10)
    plt.ylabel('arcsec', size=10)
    plt.title('noise', size=12)
    plt.savefig(plotname)

    os.system('cp ' + plotname + ' ' + public_dir + '/.')
    os.system('cp ' + data_fits + ' ' + public_dir + '/.')

    return

#_____________________________

def plot_uv(kw,current_time, ska_par):
    # make plots of the uv data
    outdir = '/tmp/ska_calculator/data/runs/' + current_time + '/'
    data = outdir + 'data.uv'
    plotname = outdir + 'uu_vv.png'

    #public_dir = 'ska_calculator/public/data/' + current_time + '/'
    public_dir = '/home/ec2-user/turbogears/ska-calculator/ska_calculator/public/data/' + current_time +'/'

    # determine the increment based on teh number of uv-points (maximum points = 16777215)
    
    if kw['tel'] == 'SKA1_Low':
        nr_baselines  = float(ska_par['low_nant']) * (float(ska_par['low_nant'])-1)/2
    if kw['tel'] == 'SKA1_Mid':
        nr_baselines  = float(ska_par['mid_nant']) * (float(ska_par['mid_nant'])-1)/2

    if kw['mode'] == 'line':
        channels = 1
    if kw['mode'] == 'cont':
        channels = float(ska_par['cont_chan'])

    cycle_time = ska_par['sample_time']
    cycle_time = float(cycle_time) / 3600 #hours

    cycles = (float(kw['hamax']) - float(kw['hamin'])) / cycle_time
    nr_uvpoints = nr_baselines * channels * cycles
    
    if nr_uvpoints > 16777215:
        increment = ceil(nr_uvpoints / 16777215)
    else:
        increment = 1

    miriad.uvplt(vis=data, axis="uu,vv", options='nobase,notitle,nofqav', select='increment('+ str(increment) + ')', device=plotname+'/png')

    os.system('cp ' + plotname + ' ' + public_dir + '/.')

    return




# ______________________________________________

def ska_sens(kw, current_time):
    # calculate ska sensitivities
    output = kw
    if kw['line'] == 'continuum':
        rest_freq = float(kw['freq'])
    if kw['line'] == 'HI':
        rest_freq = 1420.405752
    if kw['line'] == 'OH1':
        rest_freq = 1612.231
    if kw['line'] == 'OH2':
        rest_freq = 1665.402
    if kw['line'] == 'OH3':
        rest_freq = 1667.359
    if kw['line'] == 'OH4':
        rest_freq = 1720.530




    int_time = float(kw['hamax']) - float(kw['hamin'])
    output['total_time'] = kw['total_time']
    output['time'] = current_time
    output['psf_png']='data/' + kw['time'] + '/psf.png'
    output['psf_fits']='data/' + kw['time'] + '/psf.fits'
    output['noise_png']='data/' + kw['time'] + '/noise.png'
    output['noise_fits']='data/' + kw['time'] + '/noise.fits'
    output['uu_vv_png']='data/' + kw['time'] + '/uu_vv.png'
    
    outdir = '/tmp/ska_calculator/data/runs/' + current_time + '/'
    psf = outdir + 'psf.im'
    noise = outdir + 'noise.im'


    # Get the beam parameters
    #pixels = 2048
    #pixels = 639
    pixels = float(miriad.gethd(In=psf + '/naxis2'))

    min_b = int(np.floor(pixels/2 - 10))
    max_b = int(np.ceil(pixels/2 + 10))
    beam_area = 'boxes(' + str(min_b) + ',' + str(min_b) + ',' + str(max_b) + ',' + str(max_b)  + ')'
    beam = miriad.imfit(In=psf, object='beam', region=beam_area)

    # beam is a large string from which the actual parameters have to be extracted
    find_beam = 0
    for i in range(len(beam)):
        if find_beam == 0:
            if beam[i:i+20] == 'Major axis (arcsec):':
                i = i + 20
                bmaj = float(beam[i:i+16])
            if beam[i:i+20] == 'Minor axis (arcsec):':
                i = i + 20
                bmin = float(beam[i:i+16])
            if beam[i:i+25] == 'Position angle (degrees):':
                i = i + 25
                bpa = float(beam[i:i+12])
                find_beam = 1

    output['bmaj'] = bmaj
    output['bmin'] = bmin
    output['bpa'] = bpa

    # Get the noise estimates
    # estimate the noise from the header
    header = miriad.prthd(In=noise)
    find_rms = 0
    for i in range(len(header)):
        if find_rms == 0:
            if header[i:i+24] == 'Nominal Theoretical Rms:':
                i = i+24
                rms = float(header[i:i+14]) * 1000 #[mJy]
                find_rms = 1
                
    

    # get the NHI estimates
    vel_width = float(kw['width'])/float(kw['freq']) * 300000 # width [MHz] freq [GHz]
    z = (rest_freq / float(kw['freq'])) - 1
    # !! Double check this equation with Martin
    Tb = 606 * rms / (bmaj * bmin) * (1+z)**3
    NHI = 1.823e18 * Tb * vel_width


    # scale the sensitivity values if required:
    total_time = float(kw['total_time'])
    if total_time > int_time:
        rms = rms/np.sqrt(total_time/int_time)
        NHI = NHI/np.sqrt(total_time/int_time)
        Tb = Tb/np.sqrt(total_time/int_time)

    output['rms'] = "{:.3g}".format(rms)
    output['vel_width'] = "{:.3g}".format(vel_width)
    output['z'] = "{:.3g}".format(z)
    output['NHI'] = "{:.3g}".format(NHI)
    output['Tb'] = "{:.3g}".format(Tb)


    # NHI only makes sense for HI line, so leave empty for others:
    if kw['line'] == 'HI':
        output['NHI_string'] = 'Column density NHI: ' + str(output['NHI']) + ' cm^-2 \n'
    else:
        output['NHI_string'] = ''
    
    # redshift only makes sense for spectral line
    if kw['line'] == 'continuum':
        output['z_string'] = ''
    else:
        output['z_string'] = 'At observing frequency ' + str(kw['freq']) + ' MHz, the redshift z=' + str(output['z'])
        
    return output

#______________________________________

def isfloat(value):
    # check wether a given input can be a flaot    
    try:
        float(value)
        return True
    except ValueError:
        return False

#______________________________


def check_input(kw, ska_par):

    message = ''
    check = 'Pass'
    # check whether SKA1_Mid is chosen (Low not yet implemented)
    #if kw['tel'] == 'SKA1_Low':
    #    check = 'Fail'
    #    message = 'Sorry, SKA_Low is not yet implemented'
    #    return check, message
    
    # check whether spectral line is chosen (continuum not yet implemented)
    #if kw['mode'] == 'cont':
    #    check = 'Fail'
    #    message = 'Sorry, continuum is not yet implemented'
    #    return check, message
    
    # check whether band1 is chosen for Low
    if kw['tel'] == 'SKA1_Low' and  kw['band'] != 'band1':
        check = 'Fail'
        message = 'You have selected ' + kw['band'] + ', SKA_Low only has band1'
        return check, message
    # check whether frequency box has good format
    if isfloat(kw['freq']) == False:
        check = 'Fail'
        message = 'Frequency has the wrong format'
        return check, message
    # check whether frequency has good range
    if kw['tel'] == 'SKA1_Low' and kw['band'] == 'band1':
        if float(kw['freq']) < float(ska_par['low_band1_min']) or float(kw['freq']) > float(ska_par['low_band1_max']):
            check = 'Fail'
            message = 'Your selected frequency (' + kw['freq'] + ') is outside the band limits: ' + str(ska_par['low_band1_min']) + '-' + str(ska_par['low_band1_max']) + 'GHz'
            return check, message

    if kw['tel'] == 'SKA1_Mid' and kw['band'] == 'band1':
        if float(kw['freq']) < float(ska_par['mid_band1_min']) or float(kw['freq']) > float(ska_par['mid_band1_max']):
            check = 'Fail'
            message = 'Your selected frequency (' + kw['freq'] + ') is outside the band limits: ' + str(ska_par['mid_band1_min']) + '-' + str(ska_par['mid_band1_max']) + 'GHz'
            return check, message

    if kw['tel'] == 'SKA1_Mid' and kw['band'] == 'band2':
        if float(kw['freq']) < float(ska_par['mid_band2_min']) or float(kw['freq']) > float(ska_par['mid_band2_max']):
            check = 'Fail'
            message = 'Your selected frequency (' + kw['freq'] + ') is outside the band limits: ' + str(ska_par['mid_band2_min']) + '-' + str(ska_par['mid_band2_max']) + 'GHz'
            return check, message

    if kw['tel'] == 'SKA1_Mid' and kw['band'] == 'band5':
        check = 'Fail'
        message = 'Sorry, the Tsys values for band 5 have not been inmplemented yet'
        if float(kw['freq']) < float(ska_par['mid_band5_min']) or float(kw['freq']) > float(ska_par['mid_band5_max']):
            check = 'Fail'
            message = 'Your selected frequency (' + kw['freq'] + ') is outside the band limits: ' + str(ska_par['mid_band5_min']) + '-' + str(ska_par['mid_band5_max']) + 'GHz'
            return check, message

    



    # check whether width has good format and value    
    if isfloat(kw['width']) == False:
        check = 'Fail'
        message = 'Width has the wrong format'
        return check, message   
    if float(kw['width']) > float(ska_par['bandwidth_max']):
        check = 'Fail'
        message = 'Your selected width exceeds the maximum : ' + str(ska_par['bandwidth_max']) + ' MHz'
        return check, message


    # check whether width is within the boundaries:
    if kw['tel'] == 'SKA1_Low' and kw['band'] == 'band1':
        if float(kw['freq'])-(float(kw['width'])/2)< float(ska_par['low_band1_min']) or float(kw['freq'])+(float(kw['width'])/2) > float(ska_par['low_band1_max']):
            check = 'Fail'
            message = 'Your selected width goes outside the band limits: ' + str(ska_par['low_band1_min']) + '-' + str(ska_par['low_band1_max']) + 'GHz'
            return check, message

    if kw['tel'] == 'SKA1_Mid' and kw['band'] == 'band1':
        if float(kw['freq'])-(float(kw['width'])/2) < float(ska_par['mid_band1_min']) or float(kw['freq'])+(float(kw['width'])/2) > float(ska_par['mid_band1_max']):
            check = 'Fail'
            message = 'Your selected width goes outside the band limits: ' + str(ska_par['mid_band1_min']) + '-' + str(ska_par['mid_band1_max']) + 'GHz'
            return check, message

    if kw['tel'] == 'SKA1_Mid' and kw['band'] == 'band2':
        if float(kw['freq'])-(float(kw['width'])/2) < float(ska_par['mid_band2_min']) or float(kw['freq'])+(float(kw['width'])/2) > float(ska_par['mid_band2_max']):
            check = 'Fail'
            message = 'Your selected width goes outside the band limits: ' + str(ska_par['mid_band2_min']) + '-' + str(ska_par['mid_band2_max']) + 'GHz'
            return check, message

    if kw['tel'] == 'SKA1_Mid' and kw['band'] == 'band5':
        check = 'Fail'
        message = 'Sorry, the Tsys values for band 5 have not been inmplemented yet'
        if float(kw['freq'])-(float(kw['width'])/2) < float(ska_par['mid_band5_min']) or float(kw['freq'])+(float(kw['width'])/2) > float(ska_par['mid_band5_max']):
            check = 'Fail'
            message = 'Your selected width goes outside the band limits: ' + str(ska_par['mid_band5_min']) + '-' + str(ska_par['mid_band5_max']) + 'GHz'



    # check whether the declination is valid
    if isfloat(kw['dec']) == False:
        check = 'Fail'
        message = 'Declination has the wrong format'
        return check, message
    if float(kw['dec']) > 90 or float(kw['dec']) < -90:
        check = 'Fail'
        message = 'The chosen declination is outside the valid range: -90 - 90 degrees'
        return check, message

    # check whether MFS is chosen in continuum mode
    if kw['mode'] == 'cont' and kw['line'] != 'continuum':
        check = 'Fail'
        message = 'You have selected the wrong line for continuum mode, Multi Frequency Synthesis (MFS) should have been selected'
        return check, message
    if kw['mode'] == 'line' and kw['line'] == 'continuum':
        check = 'Fail'
        message = 'You have not selected a line-frequency for spectral line mode'
        


    # check whether the integration time is valid
    if isfloat(kw['total_time']) == False:
        check = 'Fail'
        message = 'The integration time has the wrong format'
        return check, message

    # check whether the hour angle is valid
    if isfloat(kw['hamin']) == False:
        check = 'Fail'
        message = 'Minimum hour angle has the wrong format'
        return check, message
    if float(kw['hamin']) > 6 or float(kw['hamin']) < -6:
        check = 'Fail'
        message = 'The chosen minimum hour angle is outside the valid range: -6 to 6 hours'
        return check, message

    if isfloat(kw['hamax']) == False:
        check = 'Fail'
        message = 'Maximum hour angle has the wrong format'
        return check, message
    if float(kw['hamax']) > 6 or float(kw['hamin']) < -6:
        check = 'Fail'
        message = 'The chosen maximum hour angle is outside the valid range: -6 to 6 hours'
        return check, message
    if float(kw['hamax']) <= float(kw['hamin']):
        check = 'Fail'
        message = 'The maximum hour angle should be larger than the minimum hour angle'
        return check, message
    
    # check whether ha range does not exceed integration time
    ha_range = float(kw['hamax']) - float(kw['hamin'])
    if ha_range > float(kw['total_time']):
        check = 'Fail'
        message = 'The selected hour angle range exceeds the total integration time'
        return check, message

    # check whether the taper is valid.
    if isfloat(kw['taper']) == False:
        check = 'Fail'
        message = 'Taper size / beam resolution has the wrong format'
        return check, message
    # calculate the allowed range
    if kw['tel'] == 'SKA1_Mid':
        qmax = float(ska_par['mid_qmax'])
    if kw['tel'] == 'SKA1_Low':
        qmax = float(ska_par['low_qmax'])
    max_res = 1.2 * (300/float(kw['freq']))/qmax * (180/np.pi) * 3600    
    if float(kw['taper']) > 1000 or float(kw['taper']) < max_res:
        check = 'Fail'
        message = 'The chosen taper is outside the valid range at this frequency: ' + str(max_res) + ' to 1000 arcsec'
        return check, message


    else:
        # there are no problems
        return check, message


