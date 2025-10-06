def payloaddatarate(rad_planet, h_orbit, gravparam, swathwidth, pixelsize, bitsperpixel):
    v_orb = math.sqrt(gravparam / (rad_planet + h_orbit))
    v_ground = v_orb * rad_planet / (rad_planet +h_orbit)
    swath_time = h_orbit * math.tan(pixelsize) / v_ground
    pixelsperswath = swathwidth / pixelsize
    pixelrate = pixelsperswath / swath_time
    datarate = pixelrate * bitsperpixel
    return datarate

def downlinkdatarate(rad_planet, h_orbit, gravparam, payloaddatarate, dutycycle, downlinktime):
    v_orb = math.sqrt(gravparam / (rad_planet + h_orbit))
    T_orb = 2 * math.pi * (rad_planet + h_orbit) / v_orb
    

