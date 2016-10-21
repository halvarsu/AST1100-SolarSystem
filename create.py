import numpy as np
import glob
import MySolarSystem as mss
import argparse

def create_pictures(everything = True):
    '''This function:
        - Loads arrays from projection folder
        - Splits filename into folders and name
        - Changes folder name to data/pictures
        - Splits name into name and extension
        - Changes extension to png
        - Changes filename from proj* to pic*
        - Recombines filename and saves array as image
        '''
    from PIL import Image
    if everything:
        i = 0
        files = glob.glob('data/projections/proj*')
        for proj_file in files:
            img = Image.fromarray(np.load(proj_file))
            pic_name = proj_file.split('/')
            pic_name[1] = 'pictures'
            fname = pic_name[-1].split('.')
            if len(fname) > 2:
                print "Something is weird with this file:",proj_file
            fname[-1] = 'png'
            fname[0] = fname[0].replace('proj', 'pic')
            pic_name[-1] = '.'.join(fname)
            pic_file = '/'.join(pic_name)
            print "Saving pic:", pic_file, "(%d/%d)" %(i,len(files))
            img.save(pic_file)
            i+= 1
            

def create_projections(start_deg, stop_deg, N):
    system = mss.MySolarSystem(87464)
    individually = True
    all_in_one = True
    super_array = np.zeros((N,480,640,3), dtype='uint8')
    super_file_name = 'all_the_projections'

    degrees = np.linspace(start_deg,stop_deg, N+1)
    phi = degrees*2*np.pi/360.

    succ_loaded = []
    for i,phi_0 in enumerate(phi[:-1]):
        name = 'data/projections/proj{:0>4}'.format(int(degrees[i]))

        #Check if file alreay exists
        if glob.glob(name+'.npy'): 
            succ_loaded.append(degrees[i])
            if all_in_one:
                proj = np.load(name+'.npy')

        else:
            if len(succ_loaded)>1:
                print "Projections at degrees %d - %d already exists, loaded files instead"%(succ_loaded[0],succ_loaded[-1])
                print '-'*40
            elif len(succ_loaded)==1:
                print "Projection at degree %d already exists, loaded file instead"%( succ_loaded[0])
                print '-'*40
            succ_loaded = []
            print 'Calculating deg: ',degrees[i], ', phi:', phi_0
            proj = system.projection(phi_0)
            print "-> Saving file "+name+'.npy\n'+'-'*40
            np.save(name, proj)
        super_array[i] = proj

    if all_in_one:
        print 'Saving superfile, ', super_file_name + '.npy'
        np.save(super_file_name, super_array)
        
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-pic", "--pictures", 
                        help="Creates pictures from arrays", 
                        action = 'store_true', default = False)
    parser.add_argument("-start", "--start_degree", 
                        help="Where to start creating", default = '0')
    parser.add_argument("-stop", "--stop_degree", 
                        help="Where to stop creating", default = '360')
    return parser.parse_args()

if __name__ == "__main__":
    args = get_args()
    start_deg = int(args.start_degree)
    stop_deg  = int(args.stop_degree)
    N = stop_deg - start_deg 
    if args.pictures:
        create_pictures()
    else:
        create_projections(start_deg, stop_deg, N)
