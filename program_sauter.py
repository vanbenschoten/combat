from __future__ import division
from cctbx.array_family import flex
from cctbx import crystal
from cctbx import miller
from iotbx.scalepack import merge

def run(args):

  vars = get_input_dict(args)
  
  if vars['xia2'] == 'yes':
    unit_cell_params = xia2(vars['bragg'])

  else:
    unit_cell_params = [vars['cella'], vars['cellb'], vars['cellc'], vars['anglea'], vars['angleb'],vars['anglec'], vars['spacegroup']]

  files = diffuse_conversion(vars['diffuse'])

  os.system("mkdir " + vars['diffuse'] +"/processed")

  for file in files:
    frame_processing(vars['diffuse'],file)

  for file in files:
    frame_averaging(vars['diffuse'] + '/processed/proc.' + file)

  reference_frame(vars['diffuse'] + '/processed/proc.', vars['reference'])

  for file in files:
    frame_scaling(vars['diffuse'] + '/processed/proc.', file)

  genlat(vars['diffuse'],files)

  args = args_generator(vars['diffuse'], vars['index_1'], vars['index_2'], vars['cella'], vars['cellb'], vars['cellc'], vars['resolution'], vars['lattice_name'], vars['processors'])

  print args

  input = 'libtbx.python genlat_labelit_sauter.py bragg=/netapp/home/vanben/lunus_processing/trypsin/bragg index_1=set_1_1_00001.cbf index_2=set_1_1_00045.cbf index_3=set_1_1_00090.cbf  diffuse=/netapp/home/vanben/lunus_processing/trypsin/diffuse indexing.data=../diffuse/set_1_1_00001.pickle indexing.data=../diffuse/set_1_1_00045.pickle indexing.data=../diffuse/set_1_1_00090.pickle codecamp.maxcell=800 xia2=yes index_only=True analyze.image=45 diffuse.lattice.resolution=2.00 cell.a=55.37 cell.b=66.64 cell.c=82.00 inputlist.fname=../diffuse/genlat.input diffuse.lattice.fname=blga_05.hkl np=8 target_cell=55.37,82.00,66.64,90,90,90 known_setting=5 distl.minimum_signal_height=5'
    
  os.system(input)

    #"args" are the arguments normally passed to integrate.py on the command line (put into get_input_dict command)
    #diffuse_integration(args)

  pre_friedel = map_symmetry_extension('out_centered.hkl', unit_cell_params)

  friedel_hkl(pre_friedel)

  
def get_input_dict(args):
# return dictionary of keywords and their values

  input_dict={}
  for arg in args:
  #print arg   
    spl=arg.split("=")
    #print spl
    if len(spl)==2:
      input_dict[spl[0]]=spl[1]

  input_dict['resolution'] = float(input_dict['resolution'])
  input_dict['processors'] = int(input_dict['processors'])
  input_dict['cella'] = float(input_dict['cella'])
  input_dict['cellb'] = float(input_dict['cellb'])
  input_dict['cellc'] = float(input_dict['cellc'])
  input_dict['anglea'] = int(input_dict['anglea'])
  input_dict['angleb'] = int(input_dict['angleb'])
  input_dict['anglec'] = int(input_dict['anglec'])

  return input_dict

def cbf2img(file):

  from dxtbx.format.Registry import Registry
  import sys, os

  f = None
  if file.split(".")[-1].lower() == "cbf" and os.path.exists(file):
    if f is None:
      f = Registry.find(file)
    img = f(file)
    db = img.get_detectorbase()
    db.readHeader()
    db.read()
    db.show_header()
    destpath = file.rstrip(".cbf") + ".img"
    print "Writing %s as %s"%(file,destpath)

    db.debug_write(destpath)

def xia2(filepath):
#Worry about changing directories?
    #os.system("cd %s"%(filepath))
    os.system("xia2 -3dii "+filepath)
    
    fin = open("xia2-summary.dat", 'r')
    lines = fin.readlines()
    unit_cell = lines[len(lines)-2]
    space_group = lines[len(lines)-1]

    q = str(unit_cell)
    a = q.split()
    del a[0]

    b = str(space_group)

    d = b.translate(None, 'Spacegroup:')
    d = d[:-1]
    d = d[1:]
    a.append(d)

    return a


def diffuse_conversion(location):
  dirListing = os.listdir(location)
  files = []
  #Better to use dictionary for this? -> Create key for each image and have corresponding
  for file in dirListing:
      if file.endswith('.cbf'):
          print file
          new_file = location + '/' + file
	  print new_file
          cbf2img(new_file)
          files.append(file.rstrip(".cbf") + ".img")
      elif file.endswith('.img'):
          files.append(file)

  return files
        

def frame_processing(filepath,file):
#frame_processing is equivalent to proc.all
    
  p = filepath + '/'
  
  #os.system("mkdir " + p +"processed")
    
  #punchim removes pixels within a specified XY region
  #os.system("punchim " + file +" 1221 2464 1242 1318 image0.img")

  #thrshim removes pixels above and below a given threshold
  os.system("thrshim " + p + file + " 1 10000 "+ p + "image.img")

  #polarim corrects for beam polarization
  os.system("polarim "+ p + "image.img " + p + "image00.img 400 0.8")


  #normim corrects for solid-angle normalization and detector-face rotation in a diffraction image
  os.system("normim " + p +  "image00.img " + p +  "image1.img")

  #modeim removes the Bragg peaks from an image by mode filtering using a specified mask size
  os.system("modeim "+ p +  "image1.img " + p + "image2.img 15 2")

  #Move the processed images to a new folder
  os.system("cp " + p +  "image2.img " + p +  "processed/proc." + file)

  #Clean up dummy images
  os.system("rm " + p + "image.img; rm " + p +  "image1.img; rm " +p +"image2.img")
  os.system("echo "+ file + " was successfully processed!")
  return

def frame_averaging(frame):
#frame_averaging is equivalent to proc.avgim

  file_prefix = frame.rstrip('.img')
  file_rf = file_prefix + ".rf"
  file_asc = file_prefix + ".asc"
  file_s_img = file_prefix + "s.img"
  file_s_sq = file_prefix + 's.sqr.rf'
  
  #avgrim calculates the average intensity vs. radius for an input image
  os.system("avgrim "+ frame +" " + file_rf)
  os.system("echo "+ frame + " avgrim works!")

  #binasc converts the average intensity vs. radius binary file to ASCII format (???)
  os.system("binasc 2 < " + file_rf + " > " + file_asc)
  os.system("echo "+ frame + " binasc works!")

  #subrfim subtracts I(r) (in this case the average radial intensity) from the diffraction image
  os.system("subrfim " + file_rf + " " + frame +" " +file_s_img)
  os.system("echo "+ frame + " subrfim works!")
  #avsqrim calculates the average *squared* pixel value vs. radius
  os.system("avsqrim " + file_s_img  + " " + file_s_sq)
  os.system("echo "+ frame + " avsqrim works!")

  #binasc converts the average intensity *squared* vs. radius binary file as: float input -> ascii output w/ leading index
  os.system("binasc 2 < " + file_s_sq + " > " + file_prefix + "s.sqr.asc")
  os.system("echo "+ frame + " binasc works!")

  #Cleaning up dummy images
  os.system("rm " + file_s_sq+"; rm "+ file_s_img)
  os.system("echo "+ frame + " file cleanup works!")

def reference_frame(prefix,file):
#reference_frame is equivalent to proc.makeref
# Can write in new variables for proc.makeref arguments
  data = file.split('.')
  image = prefix + data[0] +'.asc'
  os.system('echo IMAGE = ' +image)
  os.system('./proc.makeref '+ image + ' 250 1000' )
  os.system("echo "+ image + " reference frame creation works!")


def frame_scaling(prefix,file):

  data = file.split('.')
  image = prefix + data[0] + '.asc'
  os.system('echo scaled image = ' + image)
  os.system('./proc.scale ' + image + ' 250 1000')
  os.system("echo "+ image + " frame scaling works!")

  return

def genlat(prefix, files):

  fout = open(prefix + '/genlat.input', 'w')

  for file in files:
    os.system('cxi.image2pickle ' + prefix + '/processed/proc.' + file)
    file_prefix = file.rstrip('.img')
    scale = open(prefix + '/processed/proc.' + file_prefix + '.scale', 'r')
    lines = scale.readlines()
    scale_new = lines[0]
    scale_newest = scale_new.rstrip('\n')
    scale.close()

    fout.write(prefix + '/processed/proc.' + file_prefix + '.pickle ' + scale_newest + '\n')

  fout.close()

  return

def args_generator(location_diffuse, index_1, index_2, cella, cellb, cellc, resolution, lattice_name, processors):

  #Will need a way to read in the known_setting parameter

  os.system('cxi.image2pickle ' + location_diffuse + '/' + index_1)
  os.system('cxi.image2pickle ' + location_diffuse + '/' + index_2)


  split_1 = index_1.split('.')
  location_1 = split_1[0]

  split_2 = index_2.split('.')
  location_2 = split_2[0]

  args = ''

  args += 'indexing.data=%s/%s.pickle ' %(location_diffuse, location_1)
  args += 'indexing.data=%s/%s.pickle ' %(location_diffuse, location_2)
  args += 'codecamp.maxcell=800 '
  args += 'index_only=True '
  args += 'analyze.image=45 '
  args += 'diffuse.lattice.resolution=%0.2f ' %resolution
  args += 'cell.a=%0.2f ' %float(cella)
  args += 'cell.b=%0.2f ' %float(cellb)
  args += 'cell.c=%0.2f ' %float(cellc)
  args += 'inputlist.fname=%s/genlat.input ' %location_diffuse
  args += 'diffuse.lattice.fname=%s ' %lattice_name
  args += 'np=%d ' %processors


  return args

def map_symmetry_extension(map, unitcell):

  from cctbx import crystal
  #Read in hkl file and populate miller array
  inf = open(map, "r")
  indices = flex.miller_index()
  i_obs = flex.double()
  sig_i = flex.double()
  for line in inf.readlines():
    assert len(line.split())==4
    line = line.strip().split()
    i_obs_ = float(line[3])#/10000 #10000 is a uniform scale factor meant to re-size all diffuse intensities (normally too large for scalepack)
    sig_i_ = math.sqrt(i_obs_) 
    #if(abs(i_obs_)>1.e-6): # perhaps you don't want zeros
    indices.append([int(line[0]),int(line[1]),int(line[2])])
    i_obs.append(i_obs_)
    sig_i.append(sig_i_)
  inf.close()

  # get miller array object
  cs = crystal.symmetry(unit_cell=(float(unitcell[0]), float(unitcell[1]), float(unitcell[2]), float(unitcell[3]), float(unitcell[4]), float(unitcell[5])), space_group=unitcell[6])
  ma = miller.array(miller_set=miller.set(cs, indices), data=i_obs, sigmas=sig_i)
  ma.set_observation_type_xray_intensity()
  ma_anom = ma.customized_copy(anomalous_flag=False)
  #print ma_anom.anomalous_flag()
  ma_p1 = ma_anom.expand_to_p1()

  scalepack.merge.write(file_name='output_sequential.sca', miller_array=ma_p1)

  return 'output_sequential.sca'

def friedel_hkl(map):

  map_new = map.rstrip('.sca')

  fin = open(map,'r')
  fout = open(map_new + '.hkl', 'w')

  lines = fin.readlines()

  for line in lines:


    data = line.split()

    if len(data) == 5:

      h = int(data[0])
      h_new = -1*h
      k = int(data[1])
      k_new = -1*k
      l = int(data[2])
      l_new = -1*l
      i = float(data[3])

      #sig = data[4]

      x_1 = (str(h)+ ' ' + str(k) + ' ' + str(l) + ' ' + str(i) + '\n')
      x_2 = (str(h_new)+ ' ' + str(k_new) + ' ' + str(l_new) + ' ' + str(i) + '\n')

      fout.write(x_1)
      fout.write(x_2)

  fin.close()
  fout.close()


if __name__=="__main__":

  import sys
  import os
  import math
  
  run(sys.argv[1:])


