from time import clock, time
import numpy as np
import re
from xfel.cxi.display_spots import run_one_index_core
from cctbx.array_family import flex
from labelit.command_line.imagefiles import QuickImage
from multiprocessing import Pool

def get_input_dict(args):
  # return dictionary of keywords and their values
  input_dict={}
  for arg in args:
    spl=arg.split("=")
    if len(spl)==2:
      input_dict[spl[0]]=spl[1]
  return input_dict

def get_image(imageno, container):
  import os
  # helper function adapts the internal state from LABELIT into an image object for study
  # not intended to be altered by user
  name_of_interest = results.organizer.Files.filenames.FN[0]
  #  temp_of_interest = os.path.join(name_of_interest.cwd,name_of_interest.template)
  temp_of_interest = os.path.join(name_of_interest.cwd,"snc_tpnm_###.img")
  cnt = temp_of_interest.count('#')
  format = "%%0%dd"%cnt
  path_of_interest = temp_of_interest.replace('#'*cnt, format%imageno)
  from labelit.command_line.imagefiles import QuickImage
  return QuickImage(path_of_interest)

def pixmap(Isize1,Isize2,this_frame_phi_deg,pixel_size,size1,size2,spot_convention,procid):
  from spotfinder.applications.xfel import cxi_phil
  from iotbx.detectors.context.spot_xy_convention import spot_xy_convention
  SXYC = spot_xy_convention(pixel_size*size1,pixel_size*size2)
  from spotfinder.math_support import pixels_to_mmPos
  chunksize = int(Isize1/nproc)
  if (Isize1 % nproc !=  0):
    chunksize += 1
  x1 = procid*chunksize
  x2 = x1 + chunksize
  if (x2>Isize1):
    x2=Isize1
    #  print "Processing ",procid,"..."
  raw_spot_input = flex.vec3_double()
  for x in xrange(x1,x2): # slow dimension
    for y in xrange(Isize2): # fast dimension
      mmPos = pixels_to_mmPos(x,y,pixel_size)
      #mmPos = pixels_to_mmPos(y,x,pixel_size)
      rawspot = (mmPos[0],mmPos[1],this_frame_phi_deg)
      transpot = SXYC.select(rawspot,spot_convention)
      raw_spot_input.append(transpot)
  return raw_spot_input

def pixmapstar(args):
  return pixmap(*args)

def procimg(Isize1,Isize2,scale,mask_tag,A_matrix,rvec,DATA,latxdim,latydim,latzdim,procid):
  #  print "In procimg() with procid = ",procid
  from scitbx.matrix import col
  i0=latxdim/2-1
  j0=latydim/2-1
  k0=latzdim/2-1
  latsize = latxdim*latydim*latzdim
  lat = np.zeros(latsize*2, dtype=np.float32).reshape((2,latsize))

  #  print "Chunking ",procid,"..."

  chunksize = int(Isize2/nproc)
  if (Isize2 % nproc !=  0):
    chunksize += 1
  y1 = procid*chunksize
  y2 = y1 + chunksize
  if (y2>Isize2):
    y2=Isize2
    #  print "Processing ",procid,"..."
  for x in xrange(Isize1): # slow dimension
    for y in xrange(y1,y2): # fast dimension
      z = x*Isize2 + y
      tmid = clock()
      #      print "Calculate H ",procid,"..."
      H = A_matrix * col(rvec[z])
      #      telmatmul += clock()-tmid
      #      print "Obtaining hkl ",procid,"..."
      if (H[0]<0):
        hh = int(H[0]-.5)
      else:
        hh = int(H[0]+.5)
      if (H[1]<0):
        kk = int(H[1]-.5)
      else:
        kk = int(H[1]+.5)
      if (H[2]<0):
        ll = int(H[2]-.5)
      else:
        ll = int(H[2]+.5)
      dh = abs(H[0]-hh)
      dk = abs(H[1]-kk)
      dl = abs(H[2]-ll)
      val = int(DATA[(x,y)])
      if ((val != mask_tag) & (val != 0) & (dh > .25)&(dk > .25)&(dl > .25)):
        i = int(H[0]+i0+.5)
        j = int(H[1]+j0+.5)
        k = int(H[2]+k0+.5)
        if ((i>0)&(j>0)&(k>0)&(i<latxdim)&(j<latydim)&(k<latzdim)):
          index = k*latxdim*latydim + j*latxdim + i
          if ((val>0) & (val < 32767)):
            lat[0][index] += val*scale
            lat[1][index] += 1
            #          else:
            #            print "image %s, slow %d, fast %d, value %d"%(imgname,x,y,val)
        #      print "slow %4d fast %4d signal %8d HKL= %6.2f %6.2f %6.2f"%(
        #        y, x, DATA[(x,y)], H[0], H[1], H[2])
        #        z+=1

  return lat

def procimgstar(args):
  #  print "In procimgstar()"
  return procimg(*args)

if __name__=="__main__":
  import sys

  args = sys.argv[1:] # normally the user puts these things on command line, not in quotes, no commas
  usage = ["indexing.data=/net/sunbird/raid1/sauter/rawdata/pilatus/ribosome/images/colD55A_13_1_00001.cbf",
          # user can input any number of indexing.data image file names
          # if more than two, maxcell (unit cell upper bound in Angstroms) must be given
          # using abutting images works but slows things down
          "indexing.data=/net/sunbird/raid1/sauter/rawdata/pilatus/ribosome/images/colD55A_13_1_00401.cbf",
          "codecamp.maxcell=800",
          "index_only=True",
          "analyze.image=201"] #image number to be used for pixel analysis.
                               # but it doesn't have to be one of the images used to index.

 # Read command line arguments

  #settingidx = [a.find("known_setting")==0 for a in args].index(True)
  #known_setting = int(args.pop(settingidx).split("=")[1])
  imageidx = [a.find("analyze.image")==0 for a in args].index(True)
  image = int(args.pop(imageidx).split("=")[1])
  nprocidx = [a.find("np")==0 for a in args].index(True)
  nproc = int(args.pop(nprocidx).split("=")[1])
  cellaidx = [a.find("cell.a")==0 for a in args].index(True)
  cella = float(args.pop(cellaidx).split("=")[1])
  cellbidx = [a.find("cell.b")==0 for a in args].index(True)
  cellb = float(args.pop(cellbidx).split("=")[1])
  cellcidx = [a.find("cell.c")==0 for a in args].index(True)
  cellc = float(args.pop(cellcidx).split("=")[1])
  residx = [a.find("diffuse.lattice.resolution")==0 for a in args].index(True)
  res = float(args.pop(residx).split("=")[1])
  ifnameidx = [a.find("inputlist.fname")==0 for a in args].index(True)
  ifname = args.pop(ifnameidx).split("=")[1]
  ofnameidx = [a.find("diffuse.lattice.fname")==0 for a in args].index(True)
  ofname = args.pop(ofnameidx).split("=")[1]
  latxdim = (int(cella/res)+1)*2
  latydim = (int(cellb/res)+1)*2
  latzdim = (int(cellc/res)+1)*2

  import os

  f = open(ifname,"r")
  lines = []
  for line in f:
    if ((line.strip()!="") & (line[0] != '.')):
      lines.append(line)
  f.close()

  from spotfinder.applications.xfel import cxi_phil
  horizons_phil = cxi_phil.cxi_versioned_extract(args)

  print "indexing..."
  t0 = clock()
  results = run_one_index_core(horizons_phil)
  tel = clock()-t0
  print "done indexing (",tel," sec)"

  latsize = latxdim*latydim*latzdim
  lat = np.zeros(latsize, dtype=np.float32)
  ct = np.zeros(latsize, dtype=np.float32)

  name_of_interest = results.organizer.Files.filenames.FN[0]
  AI = results.indexing_ai
  i0=latxdim/2-1
  j0=latydim/2-1
  k0=latzdim/2-1
  mask_tag = 32767

  fileidx = 0

  #Create multiprocessor pool

  pool = Pool(processes=nproc)

  for line in lines:

    words = line.split()

    imgname = words[0]
    scale = float(words[1])

    # path_of_interest = os.path.join(name_of_interest.cwd,imgname)

    # print path_of_interest

    print "processing file %s with scale factor %f"%(imgname,scale)
    I = QuickImage(imgname)
    #    fileidx += 1
    # I = get_image(fileidx,results)
    I.read()
    DATA = I.linearintdata
  
    print "transform pixel numbers to mm positions and rotational degrees"
    from iotbx.detectors.context.spot_xy_convention import spot_xy_convention
    SF = results.spotfinder_results
    SXYC = spot_xy_convention(SF.pixel_size*SF.size1,SF.pixel_size*SF.size2)
    from spotfinder.math_support import pixels_to_mmPos
    this_frame_phi_deg = I.deltaphi/2.0+I.osc_start

    print "Creating pixel map in parallel..."
    t0 = clock()
    raw_spot_input_all = flex.vec3_double()
    Isize1 = I.size1
    Isize2 = I.size2
    pixel_size = SF.pixel_size
    pixmap_tasks = [(Isize1,Isize2,this_frame_phi_deg,SF.pixel_size,SF.size1,SF.size2,results.horizons_phil.spot_convention,procid) for procid in range(nproc)]
    raw_spot_input_it = pool.map(pixmapstar,pixmap_tasks)
    for raw_spot_input_this in raw_spot_input_it:
      raw_spot_input_all.extend(raw_spot_input_this)
      #      for j in raw_spot_input_it[i]:
      # raw_spot_input_all.append(j)
      # print "len(raw_spot_input_all) = ",len(raw_spot_input_all),"; I.size1*I.size2 = ",I.size1*I.size2
    tel = clock()-t0
    print "done creating pixel map (",tel," sec)"
    

    print "transform to laboratory axis reciprocal space coordinates"
    AI.setData(raw_spot_input_all)
    f = AI.film_to_camera()
    rvec = AI.camera_to_xyz()
  
    print "transform to fractional miller indices and populate diffuse lattice"
    from scitbx.matrix import sqr,col
    largest_setting = max(m["counter"] for m in results.pd["lattice_characters"])
    mycharacter = [m for m in results.pd["lattice_characters"] if m["counter"]==largest_setting][0]
    A_matrix = sqr(mycharacter["orient"].direct_matrix())

    #from scitbx.matrix import sqr,col
    #mycharacter = [m for m in results.pd["lattice_characters"] if m["counter"]==known_setting][0]
    #A_matrix = sqr(mycharacter["orient"].direct_matrix())

    print "Integrating diffuse scattering in parallel..."
    telmatmul=0
    t0 = clock()
    #    z = 0
    Isize1 = I.size1
    Isize2 = I.size2
    tasks = [(Isize2,Isize2,scale,mask_tag,A_matrix,rvec,DATA,latxdim,latydim,latzdim,procid) for procid in range(nproc)]
    latit = pool.map(procimgstar,tasks)
    tel = clock()-t0
    print "done integrating diffuse scattering (",tel," sec wall clock time)"
    t0 = clock()
    for l in latit:
      lat = np.add(lat,l[0])
      ct = np.add(ct,l[1])
    tel = clock()-t0
    print "Took ",tel," secs to update the lattice"
    
  for index in range(0,latsize):
    if ((ct[index] > 0) & (lat[index] != mask_tag)):
      lat[index] /= ct[index]
    else:
      lat[index] = 0


  vtkfile = open(ofname + '_raw.vtk',"w")

  a_recip = 1./cella
  b_recip = 1./cellb
  c_recip = 1./cellc


  print >>vtkfile,"# vtk DataFile Version 2.0"
  print >>vtkfile,"Generated using labelit tools"
  print >>vtkfile,"ASCII"
  print >>vtkfile,"DATASET STRUCTURED_POINTS"
  print >>vtkfile,"DIMENSIONS %d %d %d"%(latxdim,latydim,latzdim)
  print >>vtkfile,"SPACING %f %f %f"%(a_recip,b_recip,c_recip)
  print >>vtkfile,"ORIGIN %f %f %f"%(-i0*a_recip,-j0*b_recip,-k0*c_recip)
  print >>vtkfile,"POINT_DATA %d"%(latsize)
  print >>vtkfile,"SCALARS volume_scalars float 1"
  print >>vtkfile,"LOOKUP_TABLE default\n"

  index = 0
  for k in range(0,latxdim):
    for j in range(0,latydim):
      for i in range(0,latzdim):
        print >>vtkfile,lat[index],
        index += 1
      print >>vtkfile,""