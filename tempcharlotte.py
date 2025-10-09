def payloaddatarate(rad_planet, h_orbit, gravparam, swathwidth, pixelsize, bitsperpixel):
    v_orb = math.sqrt(gravparam / (rad_planet + h_orbit))
    v_ground = v_orb * rad_planet / (rad_planet +h_orbit)
    swath_time = h_orbit * math.tan(pixelsize) / v_ground
    pixelsperswath = swathwidth / pixelsize
    pixelrate = pixelsperswath / swath_time
    datarate = pixelrate * bitsperpixel
    return datarate

def downlinkdatarate(payloaddatarate, dutycycle, downlinktime):
    downlinkdatarate = 24 * 60 * 60 * dutycycle / 100 * payloaddatarate / (downlinktime * 60 * 60)
    return downlinkdatarate
    

