import numpy as np
def read_swc(input_name, output_name):
    """
    Function that reads swc files
    
    Args:
        neuron_file_name (str) : name of the swc file
    
    Returns:
        soma     (list) : [x y z r]
        segments (list) : [ [[px py pz pr], [cx cy cz cr]],  [[...], [...]] ] 
                          segment begins at parent's node x y z r and 
                          ends at child's node x y z r
    """

    in_  = open(input_name, 'r')
    out_ = open(output_name, 'w')
    lines = in_.readlines()
    point_dict = {}
    segments = []
    num_total_segments = 0
    num_segs_limit = 0
    soma = []
    for l in lines:
        l_ = l.strip()
        if len(l_) > 5:
            if l_[0] != "#":
                fields = l_.split()
                point_dict[fields[0]] = fields
                # Test that the first index is 1
                if fields[0] == '1' and fields[-1] == '-1':
                    soma       = fields[2:]
                    soma[3]    = 14.5e-3 
                    fields[2:] = soma
                    fields = [str(f) for f in fields]
                    point_dict[fields[0]] = fields
                    out_.write(' '.join(fields) + '\n')
        else:
            out_.write(l)
    point_keys = point_dict.keys()
    # Next create the list of segments - one for each child that has a parent
    for k in point_keys:
        child_fields = point_dict[k]
        # The child_fields[6] is the idx to which the point is linked
        if child_fields[6] in point_keys:
            if int(child_fields[-1]) == 1 :
                child = np.array([float(c) for c in child_fields[2:5]])
                soma  = np.array([float(s) for s in soma[:3]])
                dir   = (child - soma) / np.linalg.norm(child - soma)
                child = soma + dir * 14.5e-3
                point_dict[k][2:5] = child.astype('str')
                point_dict[k][5]   = str(0.57e-3)
                out_.write(' '.join(point_dict[k]) + '\n')
            # Don't take into account the segment that are connected to soma
            # Indeed, there is a intermediate point at the soma surface, for esthetic reasons
            elif not (int(child_fields[-1])==1 or int(child_fields[-1])==-1):
                parent_fields = point_dict[child_fields[6]]
                # Parent point x,y,z and r=radius
                px = float(parent_fields[2]) 
                py = float(parent_fields[3])
                pz = float(parent_fields[4]) 
                pr = 0.57e-3  
                parent = np.array([px, py, pz])
                # Child point x,y,z and r=radius
                cx = float(child_fields[2]) 
                cy = float(child_fields[3]) 
                cz = float(child_fields[4]) 
                child = np.array([cx, cy, cz])
                cr = 0.57e-3 

                dir = (child - parent) / np.linalg.norm(child - parent)
                child = parent + dir * 79.8e-3
                cx    = child[0]
                cy    = child[1]
                cz    = child[2]
                point_dict[k][2:5] = child.astype('str')
                point_dict[k][5]   = str(0.57e-3)
                out_.write(' '.join(point_dict[k]) + '\n')



input_name  = '/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/instructions/ISMRM24/n5.swc'
output_name = '/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/instructions/ISMRM24/n5_mod.swc'
read_swc(input_name, output_name)