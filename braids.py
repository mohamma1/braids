import maya.mel as mel
import maya.cmds as cmds # For maya commands
import string
import pdb
import math
import sys

DEBUG=0
FLOAT_DIFF_TOLERANCE = 0.0001 # If the difference between two floats is less than FLOAT_DIFF_TOLERANCE, we consider them equal


BDNA_PITCH=3.32 # Length in nanometers of one turn of B-DNA
BDNA_RISE=0.332 
BDNA_BASES_PER_TURN=10.5
BDNA_RADIUS = 1
FACIAL_ORIENTATION = 0 # 0 = clockwise, 1 = counter-clockwise
#SEGMENTATION_DIVISOR_MIN = 10 
#SEGMENTATION_DIVISOR_MAX = 150

#segmentation_divisor = 80
uniquimer_export_filename = ""

def debugprint(*args):
    if DEBUG == 1:
        print(args)

def norm(vec): 
    #Excpetions: if v is zero vector?
    ret = 0 
    for i in range(0,len(vec)):
        ret += math.pow(vec[i],2)
    ret = math.sqrt(ret)
    return ret

def unit_vec(vec): # unit vector along vec
    #Exceptions: if v is zero vector?
    ret = list(vec)
    denom = norm(vec)
    for i in range(0,len(vec)):
        ret[i] /= denom
        
    return ret
    
def diff(v, u): # v - u
    #Exceptions: if len of u != len of v
    ret = list(u)
    for i in range(0, len(u)):
        ret[i] = v[i] - u[i]
    return ret  

def sum(u, v): # v - u
    #Exceptions: if len of u != len of v
    ret = list(u)
    for i in range(0, len(u)):
        ret[i] = u[i] + v[i]
    return ret

def scal_mult(a, v):
    ret = list(v)
    for i in range(0, len(v)):
        ret[i] *= a
    return ret

def cross3(u,v): # cross product of two three-dimensional vectors
    #Exceptions if len(u) != 3 or len of v != 3
    ret = list(u)
    ret[0] = u[1]*v[2] - u[2]*v[1]
    ret[1] = u[2]*v[0] - u[0]*v[2]
    ret[2] = u[0]*v[1] - u[1]*v[0]
    return ret

def dot(u,v): 
    #Exceptions: if len of u != len of v
    ret = 0
    for i in range(0, len(u)):
        ret += u[i]*v[i]
    return ret
    
def ang(u,v):
    #Exceptions: len of u != len of v, u or v are zero vecs
    ret = dot(u,v)
    ret /= norm(u)
    ret /= norm(v)
    if ret - 1 > 0 and ret - 1 < FLOAT_DIFF_TOLERANCE: # when the angle is zero but cos of angle between u and v is slightly above one due to floating point accuracy issues, return 0
        return 0
    elif  ret + 1 < 0 and -1 * (ret + 1) < FLOAT_DIFF_TOLERANCE: # Like above but for angle 180
        return math.pi
    ret = math.acos(ret)
    return ret
        
def incident_edges(vtx_str):
    cmds.select(vtx_str)
    debugprint("Edges connected to vertex " + vtx_str)
    polyInfo_ret_strs = cmds.polyInfo(vertexToEdge=True) # TODO: Does polyInfo always list edges incident to a vertex in the counterclockwise order?
    #print(polyInfo_ret_strs[0])
    semicol_loc = string.find(polyInfo_ret_strs[0], ":")
    incident_edges_str = polyInfo_ret_strs[0][semicol_loc+1:len(polyInfo_ret_strs[0])]
    #print(incident_edges_str)
    incident_edges = incident_edges_str.split()
    debugprint(incident_edges)
    return incident_edges

def endpoints(edge_str):
    cmds.select(edge_str)
    debugprint("The endpoints of the edge " + edge_str)
    polyInfo_ret_strs = cmds.polyInfo(edgeToVertex=True)
    #print(polyInfo_ret_strs[0])
    semicol_loc = string.find(polyInfo_ret_strs[0], ":")
    endpts_str = polyInfo_ret_strs[0][semicol_loc+1:len(polyInfo_ret_strs[0])]
    #print(endpts_str)
    endpts_ids = endpts_str.split()
    debugprint(endpts_ids)
    return endpts_ids

def offset(vtx_str):
    mesh = vtx_str[0:string.find(vtx_str, ".")]
    vtx_posvec = cmds.xform(vtx_str, query=True, worldSpace=True, translation=True)
    vtx_incident_edges = list(incident_edges(vtx_str))
    other_endpts_pos = []
    avrg_ang = 0;
    for indx in range(0, len(vtx_incident_edges)):
        eid = vtx_incident_edges[indx]
        edge_str = mesh + ".e[" + eid + "]"
        endpts_ids = list(endpoints(edge_str))
        if (mesh + ".vtx[" + str(endpts_ids[0]) + "]" == vtx_str):
            other_endpt_str = mesh + ".vtx[" + str(endpts_ids[1]) + "]"
        else:
            other_endpt_str = mesh + ".vtx[" + str(endpts_ids[0]) + "]"
        other_endpt_posvec = cmds.xform(other_endpt_str, query=True, worldSpace=True, translation=True)
        
        next_eid = vtx_incident_edges[(indx + 1) % len(vtx_incident_edges)]
        next_edge_str = mesh + ".e[" + next_eid + "]"
        next_edge_endpts_ids = list(endpoints(next_edge_str))
        if (mesh + ".vtx[" + str(next_edge_endpts_ids[0]) + "]" == vtx_str):
            next_edge_other_endpt_str = mesh + ".vtx[" + str(next_edge_endpts_ids[1]) + "]"
        else:
            next_edge_other_endpt_str = mesh + ".vtx[" + str(next_edge_endpts_ids[0]) + "]"        
        next_edge_other_endpt_posvec = cmds.xform(next_edge_other_endpt_str, query=True, worldSpace=True, translation=True)
        avrg_ang += ang(diff(other_endpt_posvec, vtx_posvec), diff(next_edge_other_endpt_posvec, vtx_posvec))
        
    avrg_ang /= len(vtx_incident_edges)
    
    return BDNA_RADIUS / math.tan(avrg_ang / 2)

def next_base(base):
    nbase = "" # return empty string if three-prime end
    source_bases = cmds.listConnections(base, type="HelixBase", source=True, destination=False)
    #print("source_bases: ", source_bases)
    source_plugs = cmds.listConnections(base, type="HelixBase", source=True, destination=False, plugs=True)
    #print("source_plugs: ", source_plugs)
    i = 0
    if source_plugs:
        for plug in source_plugs:
            #print("plug: ", plug)
            if plug.find("backward") != -1: 
                nbase = source_bases[i]
                break
            i += 1
    #print("nbase: ", nbase)
    return nbase


def paired_base(base):
    pbase = "" # return empty string if unpaired
    connected_bases = cmds.listConnections(base, type="HelixBase")
    #print("connected_bases: ", connected_bases)
    plugs = cmds.listConnections(base, type="HelixBase", plugs=True)
    #print("plugs: ", plugs)
    i = 0
    if connected_bases:
        for plug in plugs:
            #print("plug: ", plug)
            if plug.find("label") != -1: 
                pbase = connected_bases[i]
                break
            i += 1
    #print("nbase: ", nbase)
    return pbase

def basetype(base):
    type="?" # unassigned by default
    '''
    dest_conns = cmds.listConnections(base, type="HelixBase", source=False, destination=True) # Apparently vhelix takes the type of the paired base if the base is a destination of the label plug
    is_dest = False
    if dest_conns:
        for dest_conn in dest_conns:
            if dest_conn.find("label") != -1:
                is_dest = True
                break
    '''
    
    is_dest = cmds.connectionInfo(base + ".label", isDestination=True)
    typeind = cmds.getAttr(base + ".label")
    if typeind == 4:
        return type
    
    if is_dest == False:
        if typeind == 0:
            type = "A"
        elif typeind == 1:
            type = "T"
        elif typeind == 2:
            type = "G"
        elif typeind == 3:
            type = "C"
    else: # take the type the Watson-Crick complement
        if typeind == 0:
            type = "T"
        elif typeind == 1:
            type = "A"
        elif typeind == 2:
            type = "C"
        elif typeind == 3:
            type = "G"

    return type    

def parent_helix(base):
    parents = cmds.listRelatives(base, parent=True)
    return parents[0] 

# Exports a vhelix scene to uniquimer xml
def export_uniquimerxml(filepath, scale=7):
    
    mel.eval("findFivePrimeEnds") # Leads to export of all strands in the scene
    fiveprimeends = cmds.ls(selection=True)
    f = open(filepath, "w")
    f.write(r"""<?xml version="1.0" encoding="UTF-8"?>""" + "\n")
    f.write(r"<root>" + "\n")
    f.write(r"""<motif type="3" scale=""" + "\"" + str(scale) + "\"" + r""">""" + "\n")
    f.write(r"<unit>" + "\n")
    f.write(r"<strands>" + "\n")
    
    # Write strands block
    print("Writing strand lengths ...")
    for fiveprimeend in fiveprimeends:
        sl = mel.eval("strandLength -base " + "\"" + fiveprimeend + "\"")
        f.write(r"""<strand length=""" + "\"" + str(sl[0]) + "\"" + r""">""" + "\n")
        f.write(r"""<bulges/>""" + "\n")
        f.write(r"</strand>" + "\n")
    f.write(r"</strands>" + "\n")
    
    # Assign the uniquimer strand id and base location for each base in the scene. This is needed for strand pairing information in Uniquimer.
    # Also compute the fixedsegements along the way
    uniquimer_strand_id = 1
    base_loc = {}
    fixedsegments = []
    num_bases = 0
    print("Assigning uniquimer ids for strands and bases and compute the fixed segments ...") 
    for fiveprimeend in fiveprimeends:
        num_bases += 1
        fixedsegment_start = -1
        fixedsegment_end = -1
        fixedsegment_seq = ""
        print("Strand " + str(uniquimer_strand_id))
        print("Strand " + str(uniquimer_strand_id) + " five-prime end: " + fiveprimeend)
        uniquimer_base_id = 1
        base_loc[fiveprimeend] = (uniquimer_strand_id,uniquimer_base_id) # five prime end is base 1 of the current strand in Uniquimer identification scheme
        
        if basetype(fiveprimeend) != "?": 
            fixedsegment_start = 1
            fixedsegment_end = 1
            fixedsegment_seq = basetype(fiveprimeend)
        
        #base_type = fiveprimeend[(string.find(fiveprimeend, "|") + 1):]
        nbase = next_base(fiveprimeend)            
        while nbase:
            num_bases += 1
            #print(nbase)
            uniquimer_base_id += 1
            base_loc[nbase] = (uniquimer_strand_id, uniquimer_base_id)
            
            
            if basetype(nbase) != "?":
                if fixedsegment_start != -1:
                    fixedsegment_end = uniquimer_base_id
                    fixedsegment_seq += basetype(nbase)
                else:                     
                    fixedsegment_start = uniquimer_base_id
                    fixedsegment_end = fixedsegment_start
                    fixedsegment_seq += basetype(nbase)
            if basetype(nbase) == "?":
                if fixedsegment_start != -1:
                    print("Fixed segment in strand " + str(uniquimer_strand_id) + " (" + str(fixedsegment_start) + "--" + str(fixedsegment_end) + "): " + fixedsegment_seq)
                    fixedsegments.append( (uniquimer_strand_id, fixedsegment_start, fixedsegment_end, fixedsegment_seq))
                    fixedsegment_start = -1
                    fixedsegment_end = -1
                    fixedsegment_seq = ""
            
            nbase = next_base(nbase)
            
        if fixedsegment_start != -1:
            print("Fixed segment in strand " + str(uniquimer_strand_id) + " (" + str(fixedsegment_start) + "--" + str(fixedsegment_end) + "): " + fixedsegment_seq)
            fixedsegments.append( (uniquimer_strand_id, fixedsegment_start, fixedsegment_end, fixedsegment_seq) )                       
        
        uniquimer_strand_id += 1
    
    
    # Write strand pairing information to strandpairs block of Uniquimer xml
    # TODO: bulges, double entry    
    f.write(r"<strandpairs>" + "\n")
    uniquimer_strand_id = 1    
    domain_added = {} # Dictionary to avoid double entry of double helical domains. 
    print("Writing strand pairing information ...")
    for fiveprimeend in fiveprimeends:
        print("Strand: " + str(uniquimer_strand_id))             
        pbase = paired_base(fiveprimeend)                
        if pbase != "":
            startpos1 = base_loc[fiveprimeend][1] # Indicates start of double helical domain
            startpos2 = base_loc[pbase][1]                        
            paired_strand = base_loc[pbase][0]
        else:
            startpos1 = -1 
            startpos2 = -1
            paired_strand = -1            
        endpos1 = startpos1
        endpos2 = startpos2
            
        nbase = next_base(fiveprimeend)
                
        while nbase:
            pbase = paired_base(nbase)
            
            if pbase != "":                
                if startpos1 != -1: # nbase is not the first base in the double helical domain                    
                    if (base_loc[pbase][0] == paired_strand): # keep updating endpos until finding end of double helical domain, TODO: bulges might cause problem here
                        endpos1 = base_loc[nbase][1]
                        endpos2 = base_loc[pbase][1]
                    else: # Doouble helical domain terminated by a paired base, i.e. start of a new domain
                        # First write strand pairing of the terminated domain                       
                        #<strandpair strand1="uniquimer_strand_id" strand2="paired_strand" startpos1="startpos1" startpos2="startpos2" endpos1="endpos1" endpos2="endpos2"/>                        
                        smaller = min(uniquimer_strand_id, paired_strand)
                        if smaller == uniquimer_strand_id:
                            smallers_startpos = startpos1
                        else:
                            smallers_startpos = endpos2                            
                        if (smaller, smallers_startpos) not in domain_added:
                            f.write(r"""<strandpair strand1=""" + "\"" + str(uniquimer_strand_id) + "\" " + r"""strand2=""" + "\"" + str(paired_strand) + "\" " + r"""startpos1=""" + "\"" + str(startpos1) + "\" " + r"""startpos2=""" + "\"" + str(startpos2) + "\" " + r"""endpos1=""" + "\"" + str(endpos1) + "\" " r"""endpos2=""" + "\"" + str(endpos2) + "\" " + """/>""" + "\n")
                            domain_added[(smaller, smallers_startpos)] = True
                            #print(domain_added)
                        # Update domain indicator variables
                        paired_strand = base_loc[pbase][0]
                        startpos1 = base_loc[nbase][1]
                        startpos2 = base_loc[pbase][1]
                        endpos1 = startpos1
                        endpos2 = startpos2                        
                else:  # nbase is the first base in the double helical domain                  
                    paired_strand = base_loc[pbase][0] # new strand
                    startpos1 = base_loc[nbase][1] # indicates start of double helical domain
                    startpos2 = base_loc[pbase][1]
                    endpos1 = startpos1
                    endpos2 = startpos2            
            else:                
                if startpos1 != -1: # Double helical domain terminated by an unpaired strand
                    # First write strand pairing of the terminated domain                       
                    #<strandpair strand1="uniquimer_strand_id" strand2="paired_strand" startpos1="startpos1" startpos2="startpos2" endpos1="endpos1" endpos2="endpos2"/>                        
                    smaller = min(uniquimer_strand_id, paired_strand)
                    if smaller == uniquimer_strand_id:
                        smallers_startpos = startpos1
                    else:
                        smallers_startpos = endpos2                            
                    if (smaller, smallers_startpos) not in domain_added:                    
                        f.write(r"""<strandpair strand1=""" + "\"" + str(uniquimer_strand_id) + "\" " + r"""strand2=""" + "\"" + str(paired_strand) + "\" " + r"""startpos1=""" + "\"" + str(startpos1) + "\" " + r"""startpos2=""" + "\"" + str(startpos2) + "\" " + r"""endpos1=""" + "\"" + str(endpos1) + "\" " r"""endpos2=""" + "\"" + str(endpos2) + "\" " + """/>""" + "\n")
                        domain_added[(smaller, smallers_startpos)] = True
                        #print(domain_added)
                    # reset domain indicator variables
                    paired_strand = -1 
                    startpos1 = -1 
                    startpos2 = -1 
                    endpos1 = -1
                    endpos2 = -1
                else:
                    pass                                                                                                  
            nbase = next_base(nbase)
        
        if startpos1 != -1: # Strand terminated by a paired nucleotide, write out the double helical domain contianing the nucleotide
            # First write strand pairing of the terminated domain                       
            #<strandpair strand1="uniquimer_strand_id" strand2="paired_strand" startpos1="startpos1" startpos2="startpos2" endpos1="endpos1" endpos2="endpos2"/>                        
            smaller = min(uniquimer_strand_id, paired_strand)
            if smaller == uniquimer_strand_id:
                smallers_startpos = startpos1
            else:
                smallers_startpos = endpos2                            
            if (smaller, smallers_startpos) not in domain_added:            
                f.write(r"""<strandpair strand1=""" + "\"" + str(uniquimer_strand_id) + "\" " + r"""strand2=""" + "\"" + str(paired_strand) + "\" " + r"""startpos1=""" + "\"" + str(startpos1) + "\" " + r"""startpos2=""" + "\"" + str(startpos2) + "\" " + r"""endpos1=""" + "\"" + str(endpos1) + "\" " r"""endpos2=""" + "\"" + str(endpos2) + "\" " + """/>""" + "\n")
                domain_added[(smaller, smallers_startpos)] = True
                #print(domain_added)
            
        uniquimer_strand_id += 1        
    f.write(r"</strandpairs>" + "\n")
    
    print("Writing fixed segments ...")
    if fixedsegments:
        f.write(r"<fixedSegments>" + "\n")
        for fixedsegment in fixedsegments:
            #<segment unitid="1" strandid="1" spos="1" epos="5" check="false">TTTTT</segment>
            f.write(r"""<segment unitid=""" + "\"1\" " + r"""strandid=""" + "\"" + str(fixedsegment[0]) + "\" " + r"""spos=""" + "\"" + str(fixedsegment[1]) + "\" " + r"""epos=""" + "\"" + str(fixedsegment[2]) + "\" " + r"""check="""+"\"false\""+""">"""+fixedsegment[3]+r"""</segment>"""+"\n")   
        f.write(r"</fixedSegments>" + "\n")
    #print("num_bases: " + str(num_bases))
    f.write(r"<maxSameLength>" + str(int(math.ceil(math.log(num_bases, 4)))) + r"</maxSameLength>")    
    f.write(r"</unit>" + "\n")
    f.write(r"</motif>" + "\n")
    f.write(r"</root>" + "\n")
    print("Done writing!")
    f.close()
    return



def fixedspacer_export_uniquimerxml(filepath, scale=7, spacer_length=3, spacer_type='T'):
    
    mel.eval("findFivePrimeEnds") # Leads to export of all strands in the scene
    fiveprimeends = cmds.ls(selection=True)
    f = open(filepath, "w")
    f.write(r"""<?xml version="1.0" encoding="UTF-8"?>""" + "\n")
    f.write(r"<root>" + "\n")
    f.write(r"""<motif type="3" scale=""" + "\"" + str(scale) + "\"" + r""">""" + "\n")
    f.write(r"<unit>" + "\n")
    f.write(r"<strands>" + "\n")
    strand_lengths = []
    
    # Find spots that need spacers
    # Assign the uniquimer strand id and base location for each base in the scene. This is needed for strand pairing information in Uniquimer.
    # Also compute the linkersegements along the way
    uniquimer_strand_id = 1
    base_loc = {}
    linkersegments = []
    print("Assigning uniquimer ids for strands and bases ...")
    for fiveprimeend in fiveprimeends:
        #print("Strand " + str(uniquimer_strand_id) + " five-prime end: " + fiveprimeend)
        sl = mel.eval("strandLength -base " + "\"" + fiveprimeend + "\"")
        strand_lengths.append(sl[0])
        cbase = fiveprimeend
        pbase = paired_base(cbase)
        nbase = next_base(cbase)
        if nbase:
            nbasep = paired_base(nbase)
        
        uniquimer_base_id = 1
        base_loc[fiveprimeend] = (uniquimer_strand_id,uniquimer_base_id) # five prime end is base 1 of the current strand in Uniquimer identification scheme
        #print("Base id: " + str(uniquimer_base_id)) 
        
        while nbase:
            #print("Current: " + cbase + ", Pair: " + pbase + ", Next: " + nbase + ", Next's pair: " + nbasep)            
            if pbase and nbasep and parent_helix(cbase) != parent_helix(nbase):
                strand_lengths[-1] += spacer_length  
                print("Tight crossover location at " + cbase)
                linkersegments.append( (uniquimer_strand_id, uniquimer_base_id + 1, uniquimer_base_id + spacer_length, spacer_type*spacer_length))
                uniquimer_base_id += (spacer_length + 1)
            else:
                uniquimer_base_id += 1
                
            base_loc[nbase] = (uniquimer_strand_id,uniquimer_base_id)
            #print("Base id: " + str(uniquimer_base_id)) 
            cbase = nbase
            pbase = paired_base(cbase)
            nbase = next_base(cbase)
            if nbase:
                nbasep = paired_base(nbase)
            
        uniquimer_strand_id += 1    
    
    # Write strands block
    print("Writing strand lengths ...")
    num_bases = 0
    for strand_length in strand_lengths:
        num_bases += strand_length
        #sl = mel.eval("strandLength -base " + "\"" + fiveprimeend + "\"")
        f.write(r"""<strand length=""" + "\"" + str(strand_length) + "\"" + r""">""" + "\n")
        f.write(r"""<bulges/>""" + "\n")
        f.write(r"</strand>" + "\n")
    f.write(r"</strands>" + "\n")
    
    
    
    # Write strand pairing information to strandpairs block of Uniquimer xml
    # TODO: bulges, double entry    
    f.write(r"<strandpairs>" + "\n")
    uniquimer_strand_id = 1    
    domain_added = {} # Dictionary to avoid double entry of double helical domains. 
    print("Writing strand pairing information ...")
    for fiveprimeend in fiveprimeends:
        #print("Strand: " + str(uniquimer_strand_id))
        #print(fiveprimeend)                
        pbase = paired_base(fiveprimeend)                
        if pbase != "":
            startpos1 = base_loc[fiveprimeend][1] # Indicates start of double helical domain
            startpos2 = base_loc[pbase][1]                        
            paired_strand = base_loc[pbase][0]
        else:
            startpos1 = -1 
            startpos2 = -1
            paired_strand = -1            
        endpos1 = startpos1
        endpos2 = startpos2
        cbase = fiveprimeend
        pbase = paired_base(cbase)
        nbase = next_base(cbase)
        if nbase:
            nbasep = paired_base(nbase)    
                
        while nbase:
            #print("Next base: " + nbase);
            if pbase and nbasep and parent_helix(cbase) != parent_helix(nbase): # cbase is base before tight crossover and nbase is after
                smaller = min(uniquimer_strand_id, paired_strand)
                if smaller == uniquimer_strand_id:
                    smallers_startpos = startpos1
                else:
                    smallers_startpos = endpos2                            
                
                if (smaller, smallers_startpos) not in domain_added:
                    #print(r"""Case A: <strandpair strand1=""" + "\"" + str(uniquimer_strand_id) + "\" " + r"""strand2=""" + "\"" + str(paired_strand) + "\" " + r"""startpos1=""" + "\"" + str(startpos1) + "\" " + r"""startpos2=""" + "\"" + str(startpos2) + "\" " + r"""endpos1=""" + "\"" + str(endpos1) + "\" " r"""endpos2=""" + "\"" + str(endpos2) + "\" " + """/>""" + "\n")
                    f.write(r"""<strandpair strand1=""" + "\"" + str(uniquimer_strand_id) + "\" " + r"""strand2=""" + "\"" + str(paired_strand) + "\" " + r"""startpos1=""" + "\"" + str(startpos1) + "\" " + r"""startpos2=""" + "\"" + str(startpos2) + "\" " + r"""endpos1=""" + "\"" + str(endpos1) + "\" " r"""endpos2=""" + "\"" + str(endpos2) + "\" " + """/>""" + "\n")
                    domain_added[(smaller, smallers_startpos)] = True
                    #print(domain_added)
                    # Update domain indicator variables
                paired_strand = base_loc[nbasep][0]
                startpos1 = base_loc[nbase][1]
                startpos2 = base_loc[nbasep][1]
                endpos1 = startpos1
                endpos2 = startpos2
                #print("startpos1 updated to: " + str(startpos1))
            else:           
                if nbasep != "":                
                    if startpos1 != -1: # nbase is not the first base in the double helical domain                    
                        if (base_loc[nbasep][0] == paired_strand): # keep updating endpos until finding end of double helical domain, TODO: bulges might cause problem here
                            endpos1 = base_loc[nbase][1]
                            endpos2 = base_loc[nbasep][1]
                        else: # Doouble helical domain terminated by a paired base, i.e. start of a new domain
                        # First write strand pairing of the terminated domain                       
                        #<strandpair strand1="uniquimer_strand_id" strand2="paired_strand" startpos1="startpos1" startpos2="startpos2" endpos1="endpos1" endpos2="endpos2"/>                        
                            smaller = min(uniquimer_strand_id, paired_strand)
                            if smaller == uniquimer_strand_id:
                                smallers_startpos = startpos1
                            else:
                                smallers_startpos = endpos2                            
                            if (smaller, smallers_startpos) not in domain_added:
                                #print(r"""Case B: <strandpair strand1=""" + "\"" + str(uniquimer_strand_id) + "\" " + r"""strand2=""" + "\"" + str(paired_strand) + "\" " + r"""startpos1=""" + "\"" + str(startpos1) + "\" " + r"""startpos2=""" + "\"" + str(startpos2) + "\" " + r"""endpos1=""" + "\"" + str(endpos1) + "\" " r"""endpos2=""" + "\"" + str(endpos2) + "\" " + """/>""" + "\n")
                                f.write(r"""<strandpair strand1=""" + "\"" + str(uniquimer_strand_id) + "\" " + r"""strand2=""" + "\"" + str(paired_strand) + "\" " + r"""startpos1=""" + "\"" + str(startpos1) + "\" " + r"""startpos2=""" + "\"" + str(startpos2) + "\" " + r"""endpos1=""" + "\"" + str(endpos1) + "\" " r"""endpos2=""" + "\"" + str(endpos2) + "\" " + """/>""" + "\n")
                                domain_added[(smaller, smallers_startpos)] = True
                                #print(domain_added)                          
                            # Update domain indicator variables
                            paired_strand = base_loc[nbasep][0]
                            startpos1 = base_loc[nbase][1]
                            startpos2 = base_loc[nbasep][1]                        
                            endpos1 = startpos1
                            endpos2 = endpos2                                                    
                    else:  # nbase is the first base in the double helical domain                  
                        paired_strand = base_loc[nbasep][0] # new strand
                        startpos1 = base_loc[nbase][1] # indicates start of double helical domain
                        startpos2 = base_loc[nbasep][1]
                        endpos1 = startpos1
                        endpos2 = startpos2          
                else:                
                    if startpos1 != -1: # Double helical domain terminated by an unpaired strand
                        # First write strand pairing of the terminated domain                       
                        #<strandpair strand1="uniquimer_strand_id" strand2="paired_strand" startpos1="startpos1" startpos2="startpos2" endpos1="endpos1" endpos2="endpos2"/>                        
                        smaller = min(uniquimer_strand_id, paired_strand)
                        if smaller == uniquimer_strand_id:
                            smallers_startpos = startpos1
                        else:
                            smallers_startpos = endpos2                            
                        if (smaller, smallers_startpos) not in domain_added:
                            #print(r"""Case C: <strandpair strand1=""" + "\"" + str(uniquimer_strand_id) + "\" " + r"""strand2=""" + "\"" + str(paired_strand) + "\" " + r"""startpos1=""" + "\"" + str(startpos1) + "\" " + r"""startpos2=""" + "\"" + str(startpos2) + "\" " + r"""endpos1=""" + "\"" + str(endpos1) + "\" " r"""endpos2=""" + "\"" + str(endpos2) + "\" " + """/>""" + "\n")                    
                            f.write(r"""<strandpair strand1=""" + "\"" + str(uniquimer_strand_id) + "\" " + r"""strand2=""" + "\"" + str(paired_strand) + "\" " + r"""startpos1=""" + "\"" + str(startpos1) + "\" " + r"""startpos2=""" + "\"" + str(startpos2) + "\" " + r"""endpos1=""" + "\"" + str(endpos1) + "\" " r"""endpos2=""" + "\"" + str(endpos2) + "\" " + """/>""" + "\n")
                            domain_added[(smaller, smallers_startpos)] = True
                            #print(domain_added)
                        # reset domain indicator variables
                        paired_strand = -1 
                        startpos1 = -1 
                        startpos2 = -1 
                        endpos1 = -1
                        endpos2 = -1                        
                    else:
                        pass         
             
            cbase = nbase
            pbase = paired_base(cbase)
            nbase = next_base(cbase)
            if nbase:
                nbasep = paired_base(nbase)                                                                                               
        
        if startpos1 != -1: # Strand terminated by a paired nucleotide, write out the double helical domain contianing the nucleotide
            # First write strand pairing of the terminated domain                       
            #<strandpair strand1="uniquimer_strand_id" strand2="paired_strand" startpos1="startpos1" startpos2="startpos2" endpos1="endpos1" endpos2="endpos2"/>                        
            smaller = min(uniquimer_strand_id, paired_strand)
            if smaller == uniquimer_strand_id:
                smallers_startpos = startpos1
            else:
                smallers_startpos = endpos2                            
            if (smaller, smallers_startpos) not in domain_added:
                #print(r"""Case D: <strandpair strand1=""" + "\"" + str(uniquimer_strand_id) + "\" " + r"""strand2=""" + "\"" + str(paired_strand) + "\" " + r"""startpos1=""" + "\"" + str(startpos1) + "\" " + r"""startpos2=""" + "\"" + str(startpos2) + "\" " + r"""endpos1=""" + "\"" + str(endpos1) + "\" " r"""endpos2=""" + "\"" + str(endpos2) + "\" " + """/>""" + "\n")            
                f.write(r"""<strandpair strand1=""" + "\"" + str(uniquimer_strand_id) + "\" " + r"""strand2=""" + "\"" + str(paired_strand) + "\" " + r"""startpos1=""" + "\"" + str(startpos1) + "\" " + r"""startpos2=""" + "\"" + str(startpos2) + "\" " + r"""endpos1=""" + "\"" + str(endpos1) + "\" " r"""endpos2=""" + "\"" + str(endpos2) + "\" " + """/>""" + "\n")
                domain_added[(smaller, smallers_startpos)] = True
                #print(domain_added)
            
        uniquimer_strand_id += 1
        
        
                
    f.write(r"</strandpairs>" + "\n")
    print("Writing fixed segments ...")
    if linkersegments:
        f.write(r"<fixedSegments>" + "\n")
        for fixedsegment in linkersegments:
            #<segment unitid="1" strandid="1" spos="1" epos="5" check="false">TTTTT</segment>
            f.write(r"""<segment unitid=""" + "\"1\" " + r"""strandid=""" + "\"" + str(fixedsegment[0]) + "\" " + r"""spos=""" + "\"" + str(fixedsegment[1]) + "\" " + r"""epos=""" + "\"" + str(fixedsegment[2]) + "\" " + r"""check="""+"\"false\""+""">"""+fixedsegment[3]+r"""</segment>"""+"\n")   
        f.write(r"</fixedSegments>" + "\n")
    #print("num_bases: " + str(num_bases))
    f.write(r"<maxSameLength>" + str(int(math.ceil(math.log(num_bases, 4)))) + r"</maxSameLength>")    
    f.write(r"</unit>" + "\n")
    f.write(r"</motif>" + "\n")
    f.write(r"</root>" + "\n")
    
    print("Done writing!")
    f.close()
    return



def export_to_uniquimer_dialog():
    global uniquimer_export_filename
    uniquimer_export_filename = cmds.fileDialog2(fileFilter="xml files (*.xml);; All Files (*.*)", dialogStyle=2)    
    if uniquimer_export_filename:
        print("Exporting strands to " + uniquimer_export_filename[0] + " ...")
        #export_uniquimerxml(uniquimer_export_filename[0])
        
        linker_nt = cmds.optionMenu('linker_nt_menu', query = True, value = True)
        
        if cmds.radioButton('const_linker_radio', query = True, select = True): 
            linker_length = cmds.intField('linker_length_field', query = True, value = True)               
            fixedspacer_export_uniquimerxml(uniquimer_export_filename[0], 7, linker_length, linker_nt)
        elif cmds.radioButton('var_linker_radio', query = True, select = True):
            mel.eval("fillStrandGaps")
            assign_polyn_to_ss_segs(linker_nt)
            export_uniquimerxml(uniquimer_export_filename[0])             
        

def assign_polyn_to_ss_segs(basetype="T"):
    mel.eval("findFivePrimeEnds") # Leads to export of all strands in the scene
    fiveprimeends = cmds.ls(selection=True)
    for fiveprimeend in fiveprimeends:
        #print(fiveprimeend)
        pbase = paired_base(fiveprimeend)
        if pbase == "":
            mel.eval("applySequence -target \"" + fiveprimeend + "\" -sequence \"" + basetype + "\"")
        nbase = next_base(fiveprimeend)
        while nbase:
            #print(nbase)
            pbase = paired_base(nbase)
            if pbase == "":
                mel.eval("applySequence -target \"" + nbase + "\" -sequence \"" + basetype + "\"")
            nbase = next_base(nbase)
                            
        

def facial_cycle_render():
    #debugprint("Segmentation divisor = " + str(segmentation_divisor))
    
    min_helix_bases = 19
    
    lsret = cmds.ls(selection=True, transforms=True)
    debugprint(lsret)
    if not lsret:
        print("Please select a mesh first!")
        return
    else:
        print("Mesh selected")
        
    mesh=lsret[0]
    num_edges = cmds.polyEvaluate(edge=True)
    num_vertices = cmds.polyEvaluate(vertex=True)
    num_bases = [0] * num_edges
    vtx_offsets = [0] * num_vertices
    
    print("Number of edges: " + str(num_edges))
    print("Number of vertices: " + str(num_vertices))
    
    print("Calculating offsets at the vertices ...")
    for vid in range(0, num_vertices):
        vtx_str = mesh + ".vtx[" + str(vid) + "]"
        vtx_offsets[vid] = offset(vtx_str)
    
    # For all edges in the selected mesh, create helix centered at the edge and along the edge axis
    print("---------------------------------------------------------------")
    print("Creating helices at the edges ...") 
    for eid in range(0, num_edges):
        print("Processing edge " + str(eid) + " ...")
        edge=mesh + ".e[" + str(eid) + "]"
        
        # Find the endpoints of the edge
        eid_endpts_ids = endpoints(edge)
        endpt1 = mesh + ".vtx[" + eid_endpts_ids[0] + "]"
        endpt2 = mesh + ".vtx[" + eid_endpts_ids[1] + "]"
        debugprint(endpt1)
        debugprint(endpt2)
        
        # Compute the center point and length of the edge
        endpt1_posvec = cmds.xform(endpt1, query=True, worldSpace=True, translation=True)
        endpt2_posvec = cmds.xform(endpt2, query=True, worldSpace=True, translation=True)
        endpt1_helixend = sum(endpt1_posvec, scal_mult(vtx_offsets[int(eid_endpts_ids[0])], unit_vec(diff(endpt2_posvec, endpt1_posvec))))
        endpt2_helixend = sum(endpt2_posvec, scal_mult(vtx_offsets[int(eid_endpts_ids[1])], unit_vec(diff(endpt1_posvec, endpt2_posvec))))        
        
        e_center_posvec = [0, 0, 0]
        e_center_posvec[0] = (endpt1_posvec[0] + endpt2_posvec[0]) / 2.0
        e_center_posvec[1] = (endpt1_posvec[1] + endpt2_posvec[1]) / 2.0
        e_center_posvec[2] = (endpt1_posvec[2] + endpt2_posvec[2]) / 2.0
        debugprint("Edge center", e_center_posvec) 
        e_length = math.sqrt(math.pow(endpt1_posvec[0] - endpt2_posvec[0],2) + math.pow(endpt1_posvec[1] - endpt2_posvec[1],2) + math.pow(endpt1_posvec[2] - endpt2_posvec[2],2))
        debugprint("Edge length", e_length)
        
        helix_center_posvec = [0,0,0]
        helix_center_posvec[0] = (endpt1_helixend[0] + endpt2_helixend[0]) / 2.0
        helix_center_posvec[1] = (endpt1_helixend[1] + endpt2_helixend[1]) / 2.0
        helix_center_posvec[2] = (endpt1_helixend[2] + endpt2_helixend[2]) / 2.0
        print("Helix center", helix_center_posvec)
        offset2offset_length = norm( diff(endpt1_helixend, endpt2_helixend))
        print("Offset to offset length", offset2offset_length)
        # Find the approximate number of full turns to render the edge and create a helix of that length positioned at the center of the edge
        num_bases[eid] = int((math.floor(math.floor(offset2offset_length / BDNA_PITCH)) * BDNA_BASES_PER_TURN)) #TODO: improved offsetting
        print("The approximate number of bases for an interger-number-of-full-turns helix not exceeding the offsets is " + str(num_bases[eid])) 
        
        if num_bases[eid] <= min_helix_bases:
            print("Edge " + str(eid) + " is too short! Rescale the mesh and try again.")
            sys.exit("Edge too short!")
                    
        create_helix_cmd = "createHelix -base " + str(num_bases[eid]) 
        mel.eval(create_helix_cmd)
        
        # Rotate the helix so that it aligns along the edge and so that its five prime end is placed on the FACIAL_ORIENTATION neighboring edge  
        # Compute alpha, beta and gamma for Euler rotations (wikipedia Euler angles page [accessed 21.02.2019 ]) and use them for rotateX, rotateY and rotateZ of Maya
        Z = [0,0,0] # Euler rotation Z axis = axis of edge
        Z = unit_vec(diff(endpt2_posvec, endpt1_posvec))
        debugprint("Z: ", Z)
        
        # Compute the FACIAL_ORIENTATION neighbor of the edge at endpoint 1
        endpt1_incident_edges = list(incident_edges(endpt1))
        edge_loc = endpt1_incident_edges.index(str(eid)) # location of current edge in the list of incident edges to endpoint 1
           
        if FACIAL_ORIENTATION == 0: # get clockwise neighboring edge TODO: considering removing this
            otherneigh_edge_loc = (edge_loc + 1) % len(endpt1_incident_edges)
            if edge_loc > 0:
                neigh_edge_loc = edge_loc - 1                
            else:
                neigh_edge_loc = len(endpt1_incident_edges) - 1
        else: # get counterclockwise neighboring edge
            neigh_edge_loc = (edge_loc + 1) %  len(endpt1_incident_edges)
            if edge_loc > 0:
                otherneigh_edge_loc = edge_loc - 1
            else:
                otherneigh_edge_loc = len(endpt1_incident_edges) - 1                
        
        neigh_edge = int(endpt1_incident_edges[neigh_edge_loc])
        otherneigh_edge = int(endpt1_incident_edges[otherneigh_edge_loc])
        debugprint("neigh_edge = ", neigh_edge)
        
        # Find endpt1's id to compare against the other endpoints of the neighboring edges
        endpt1_id = int(endpt1.split("[")[1].split("]")[0])
                    
        # Compute the unit vector from endpoint 1 to the other endpoint of the clockwise neighboring edge
        cmds.select(mesh + ".e[" + str(neigh_edge) + "]")
        polyInfo_ret_strs = cmds.polyInfo(edgeToVertex=True)
        debugprint(polyInfo_ret_strs)
        neigh_endpts_str = polyInfo_ret_strs[0][string.find(polyInfo_ret_strs[0], ":")+1:len(polyInfo_ret_strs[0])]
        neigh_endpts = neigh_endpts_str.split()[0:2]
        if int(neigh_endpts[0]) == endpt1_id:
            neigh_other_endpt_id = int(neigh_endpts[1])
        else:
            neigh_other_endpt_id = int(neigh_endpts[0])
        debugprint("Neighboring edges other endpoint: " + str(neigh_other_endpt_id))
        neigh_other_endpt = mesh + ".vtx[" + str(neigh_other_endpt_id) + "]"
        neigh_other_endpt_posvec = cmds.xform(neigh_other_endpt, query=True, worldSpace=True, translation=True)
        endpt1_to_neigh_other_endpt_unitvec = unit_vec(diff(neigh_other_endpt_posvec, endpt1_posvec))
        debugprint("Unit vector from " + str(endpt1_id) + " to " + str(neigh_other_endpt_id), endpt1_to_neigh_other_endpt_unitvec)
        
        # Compute the unit vector from endpoint 1 to the other endpoint of the other (counterclockwise) neighboring edge
        cmds.select(mesh + ".e[" + str(otherneigh_edge) + "]")
        polyInfo_ret_strs = cmds.polyInfo(edgeToVertex=True)
        debugprint(polyInfo_ret_strs)
        otherneigh_endpts_str = polyInfo_ret_strs[0][string.find(polyInfo_ret_strs[0], ":")+1:len(polyInfo_ret_strs[0])]
        otherneigh_endpts = otherneigh_endpts_str.split()[0:2]
        if int(otherneigh_endpts[0]) == endpt1_id:
            otherneigh_other_endpt_id = int(otherneigh_endpts[1])
        else:
            otherneigh_other_endpt_id = int(otherneigh_endpts[0])
        
        debugprint("otherneigh_edge = ", otherneigh_edge)
        debugprint("Other neighboring edges other endpoint: " + str(otherneigh_other_endpt_id))                        
        
        otherneigh_other_endpt = mesh + ".vtx[" + str(otherneigh_other_endpt_id) + "]"
        otherneigh_other_endpt_posvec = cmds.xform(otherneigh_other_endpt, query=True, worldSpace=True, translation=True)
        endpt1_to_otherneigh_other_endpt_unitvec = unit_vec(diff(otherneigh_other_endpt_posvec, endpt1_posvec))
        debugprint("Unit vector from " + str(endpt1_id) + " to " + str(otherneigh_other_endpt_id), endpt1_to_otherneigh_other_endpt_unitvec)
                
        # Compute the desired direction from the helix axis to the first forward base (5' base)
        if abs(abs(dot(Z,endpt1_to_neigh_other_endpt_unitvec)) - 1) > FLOAT_DIFF_TOLERANCE:
            debugprint("Edge " + str(eid) + " and " + str(neigh_edge) + " are neither parallel nor antiparallel" )            
            Y = list(unit_vec(cross3(Z, cross3(endpt1_to_neigh_other_endpt_unitvec, Z)))) #TODO: why is cross product of two perpendicular unit vecs not resulting in a unit vector
        else:
            debugprint("Edge " + str(eid) + " and " + str(neigh_edge) + " are parallel or antiparallel" )            
            Y = list(unit_vec(cross3(Z, cross3(scal_mult(-1,endpt1_to_otherneigh_other_endpt_unitvec), Z)))) #TODO: why is cross product of two perpendicular unit vecs not resulting in a unit vector
        
        # Compute the vertex-normal of endpt1
        cmds.select(endpt1)
        vfns = cmds.polyNormalPerVertex( query=True, xyz=True ) #vertex-face normals
        debugprint("vns: ", vfns)
        vn = [0, 0, 0] # vertex normal = average of vertex-face normals        
        for i in range(0, len(vfns)):
            if (i % 3) == 0:
                vn[0] += vfns[i]
            elif (i % 3) == 1:
                vn[1] += vfns[i]
            else:
                vn[2] += vfns[i]
        vn[0] /= (len(vfns) / 3)
        vn[1] /= (len(vfns) / 3)
        vn[2] /= (len(vfns) / 3)                 
           
        if dot(endpt1_to_neigh_other_endpt_unitvec, cross3(Z, vn)) < 0:
            Y = scal_mult(-1, Y)
        #if neigh_edge == otherneigh_edge:
            
                    
        X = unit_vec(cross3(Y, Z))
        debugprint("X, norm(X): ", X, norm(X))
        debugprint("Y, norm(Y): ", Y, norm(Y))
        debugprint("Z, norm(Z): ", Z, norm(Z))
        hid = eid + 1 # helix id, vhelix starts helix ids from 1
        cmds.select('helix' + str(hid))
        
        # Rotate and translate helix 
        cmds.xform(m = (X[0], X[1], X[2], 0, Y[0], Y[1], Y[2], 0, Z[0], Z[1], Z[2], 0, e_center_posvec[0], e_center_posvec[1], e_center_posvec[2], 1))
        #cmds.xform(scale=[1,1,1])
    #end for 
    
    # Connect helices to form the vertices
    print("---------------------------------------------------------------")
    print("Connecting the helices to form the vertices ...")
    cmds.select(mesh)
    for vid in range(0, num_vertices):
        vtx_str = mesh + ".vtx[" + str(vid) + "]"
        vtx_incident_edges = list(incident_edges(vtx_str))
        
        # TODO: Assuming maya orders the incident edges of a vertex in a manifold mesh in a counterclockwise order ...
        i = 0
        for eid in vtx_incident_edges:
            # Connect the five prime of the next edge and the three prime of current edge
            eid_endpts_ids = list(endpoints(mesh + ".e[" + eid + "]"))
            if vid == int(eid_endpts_ids[0]):
                cur_edge_3p = "helix" + str( int(eid) + 1) + "|backw_1"
            else:
                cur_edge_3p = "helix" + str( int(eid) + 1) + "|forw_" + str(num_bases[int(eid)])
    
            nxt_edge_id = int(vtx_incident_edges[(i + 1) % len(vtx_incident_edges)])
            nxt_edge_endpts_ids = list(endpoints(mesh + ".e[" + str(nxt_edge_id) + "]")) #TODO list constructor first two elements only
            
            if vid == int(nxt_edge_endpts_ids[0]):
                nxt_edge_5p = "helix" + str( nxt_edge_id + 1) + "|forw_1"
            else:
                nxt_edge_5p = "helix" + str( nxt_edge_id + 1) + "|backw_" + str(num_bases[int(nxt_edge_id)])
                    
            print("Connecting bases " + cur_edge_3p + " and " +  nxt_edge_5p + " ... ")
            connect_bases_cmd = "connectBases  -first " + cur_edge_3p + " -second " +  nxt_edge_5p
            mel.eval(connect_bases_cmd)
            i += 1
            
    # Nick the circular strands so that each of the resulting linear strands span only two edges
    print("---------------------------------------------------------------")
    print("Nicking the circular strands so that each of the resulting linear strands span only two edges ...")
    for eid in range(0, num_edges):
        
        if abs(num_bases[eid] - 21)  <= 2: # ~2 turns
            forw_disconn_base_id = num_bases[eid] - 7
            backw_disconn_base_id = 7
        elif abs(num_bases[eid] - 32)  <= 2: # ~3 turns
            forw_disconn_base_id = num_bases[eid] - 11
            backw_disconn_base_id = 11
        elif abs(num_bases[eid] - 42)  <= 2: # ~4 turns
            forw_disconn_base_id = num_bases[eid] - 14
            backw_disconn_base_id = 14
        elif num_bases[eid] > 44 and num_bases[eid] < 82: # ~ 5, 6, or 7 turns
            forw_disconn_base_id = 21
            forw_disconn_base_id2 = num_bases[eid] - 11
            backw_disconn_base_id = 11
            backw_disconn_base_id2 = num_bases[eid] - 21
        elif abs(num_bases[eid] - 84) < 2:
            forw_disconn_base_id = 21
            forw_disconn_base_id2 = 47
            forw_disconn_base_id3 = num_bases[eid] - 11
            backw_disconn_base_id = 11
            backw_disconn_base_id2 = 37
            backw_disconn_base_id3 = num_bases[eid] - 21
        elif abs(num_bases[eid] - 95) < 2:
            forw_disconn_base_id = 21
            forw_disconn_base_id2 = 52
            forw_disconn_base_id3 = num_bases[eid] - 11
            backw_disconn_base_id = 11
            backw_disconn_base_id2 = 42
            backw_disconn_base_id3 = num_bases[eid] - 21
        elif num_bases[eid] > 97:
            forw_disconn_base_id = 21
            forw_disconn_base_id2 = 21 + int((num_bases[eid] - 32) / 2) 
            forw_disconn_base_id3 = num_bases[eid] - 11
            backw_disconn_base_id = 11
            backw_disconn_base_id2 = 11 + int((num_bases[eid] - 32) / 2)
            backw_disconn_base_id3 = num_bases[eid] - 21
            
            
        forw_disconn_base = "helix" + str(eid+1) + "|forw_" + str(forw_disconn_base_id)
        backw_disconn_base = "helix" + str(eid+1) + "|backw_" + str(backw_disconn_base_id + 1)
        mel.eval("disconnectBase -target " + forw_disconn_base)
        mel.eval("disconnectBase -target " + backw_disconn_base)
        if num_bases[eid] > 44:
            forw_disconn_base2 = "helix" + str(eid+1) + "|forw_" + str(forw_disconn_base_id2)
            backw_disconn_base2 = "helix" + str(eid+1) + "|backw_" + str(backw_disconn_base_id2)
            mel.eval("disconnectBase -target " + forw_disconn_base2)
            mel.eval("disconnectBase -target " + backw_disconn_base2)
        if num_bases[eid] > 82:
            forw_disconn_base3 = "helix" + str(eid+1) + "|forw_" + str(forw_disconn_base_id3)
            backw_disconn_base3 = "helix" + str(eid+1) + "|backw_" + str(backw_disconn_base_id3)
            mel.eval("disconnectBase -target " + forw_disconn_base3)
            mel.eval("disconnectBase -target " + backw_disconn_base3)
                                            
        '''
        num_nicks = int(math.ceil(int(num_bases[eid]) / int(segmentation_divisor))) + 1
        debugprint("Number of nicks in helix " + str(eid + 1) + ": " + str(num_nicks))
        seg_size = int(num_bases[eid] / (2 * num_nicks + 1))
        for i in range(0, num_nicks):
            forw_disconn_base_id = 2 * seg_size * (i+1)
            forw_disconn_base = "helix" + str(eid+1) + "|forw_" + str(forw_disconn_base_id + 1)
            backw_disconn_base_id = seg_size + 2 * seg_size * i
            backw_disconn_base = "helix" + str(eid+1) + "|backw_" + str(backw_disconn_base_id + 1)
            forw_strand_segm_cmd = "disconnectBase -target " + forw_disconn_base
            backw_strand_segm_cmd = "disconnectBase -target " + backw_disconn_base
            print("Nicking the forward strand of helix " + str(eid + 1) + " at " + forw_disconn_base)
            print("Nicking the backward strand of helix " + str(eid + 1) + " at " + backw_disconn_base)
            mel.eval(forw_strand_segm_cmd)
            mel.eval(backw_strand_segm_cmd)
        '''
        '''
        forw_strand_segm_ratio = 2.0 / 3.0
        backw_strand_segm_ratio = 1.0 / 3.0
        debugprint("Number of bases in helix " + str(eid + 1) + " = " + str(num_bases[eid]))
        #print(str(float(num_bases[eid]) * forw_strand_segm_ratio))
        forw_disconn_base = "helix" + str(eid+1) + "|forw_" + str(int(math.floor(float(num_bases[eid]) * forw_strand_segm_ratio)))
        backw_disconn_base = "helix" + str(eid+1) + "|backw_" + str(int(math.floor(float(num_bases[eid]) * backw_strand_segm_ratio)))
        print("Nicking the forward strand of helix " + str(eid + 1) + " at " + forw_disconn_base)
        print("Nicking the backward strand of helix " + str(eid + 1) + " at " + backw_disconn_base)
        forw_disconn_base_shape = cmds.ls(forw_disconn_base, dagObjects=True, shapes=True)
        forw_disconn_base_orig_color = cmds.listConnections(forw_disconn_base_shape, type="shadingEngine")
        backw_disconn_base_shape = cmds.ls(backw_disconn_base, dagObjects=True, shapes=True)
        backw_disconn_base_orig_color = cmds.listConnections(backw_disconn_base_shape, type="shadingEngine")
        print("Pre-nicking color of strand containing base " + forw_disconn_base + ": " + forw_disconn_base_orig_color[0])
        print("Pre-nicking color of strand containing base " + backw_disconn_base + ": " + backw_disconn_base_orig_color[0])
        
        forw_strand_segm_cmd = "disconnectBase -target " + forw_disconn_base
        backw_strand_segm_cmd = "disconnectBase -target " + backw_disconn_base
        print(forw_strand_segm_cmd)
        print(backw_strand_segm_cmd)
        mel.eval(forw_strand_segm_cmd)
        mel.eval(backw_strand_segm_cmd)
        
        forw_disconn_base_new_color = cmds.listConnections(forw_disconn_base_shape, type="shadingEngine")
        backw_disconn_base_new_color = cmds.listConnections(backw_disconn_base_shape, type="shadingEngine")
        print("Post-nicking color of strand containing base " + forw_disconn_base + ": " + forw_disconn_base_new_color[0])
        print("Post-nicking color of strand containing base " + backw_disconn_base + ": " + backw_disconn_base_new_color[0])
        
        print("Repainting strand to original color ...")
        while forw_disconn_base_new_color[0] != forw_disconn_base_orig_color[0]:
            mel.eval("paintStrand -target " + forw_disconn_base)
            forw_disconn_base_new_color = cmds.listConnections(forw_disconn_base_shape, type="shadingEngine")
            #print(forw_disconn_base_new_color)
            
        while backw_disconn_base_new_color[0] != backw_disconn_base_orig_color[0]:
            mel.eval("paintStrand -target " + backw_disconn_base)
            backw_disconn_base_new_color = cmds.listConnections(backw_disconn_base_shape, type="shadingEngine")
            #print(backw_disconn_base_new_color)
        '''

            
    
    cmds.select(clear=True)
    #mel.eval("fillStrandGaps")
    #assign_polyn_to_ss_segs("T")
    
    return
    
#def set_linker_nucleotide():
#    linker_nucleotide = cmds.optionMenu(linker_nt_menu, query=True, value=True)

'''
def set_segmentation_divisor():
    global segmentation_divisor
    segmentation_divisor = cmds.intSliderGrp("segmentation_divisor_slidergrp", query=True, value=True)
    debugprint("In segmentatation divisor callback function, segmentation_divisor = " + str(segmentation_divisor))
    #if cmds.intSliderGrp("segmentation_divisor_slidergrp", exists=True):
    return
'''

win_id = "facial_cycle_renderer"
if cmds.window(win_id, exists=True):
    cmds.deleteUI(win_id)
    
cmds.window(win_id, title="Cadwire: DNA design of surface mesh wireframes " )

cmds.columnLayout(rowSpacing=10, columnWidth=200)
cmds.frameLayout( label='Helix generation', width=200)
#cmds.intSliderGrp( "segmentation_divisor_slidergrp", field=True, label="Segmentation divisor: ", minValue=SEGMENTATION_DIVISOR_MIN, maxValue=SEGMENTATION_DIVISOR_MAX, fieldMinValue=SEGMENTATION_DIVISOR_MIN, fieldMaxValue=SEGMENTATION_DIVISOR_MAX, value=80, changeCommand="set_segmentation_divisor()"  )
cmds.button( label="Generate helices", command="facial_cycle_render()")
cmds.setParent('..')
cmds.frameLayout( label='Uniquimer export', width=200 )
cmds.rowLayout(numberOfColumns=3)
cmds.text(label='Linker style: ')
linker_style_radio_coll = cmds.radioCollection()
cmds.radioButton( 'const_linker_radio', label='Constant', select = True )
cmds.radioButton( 'var_linker_radio', label='Variable' )
#cmds.radioCollection( linker_style_radio_coll, edit=True, select=const_radio )
cmds.setParent('..')
cmds.optionMenu('linker_nt_menu', label='Linker nucleotide: ' )
cmds.menuItem( label='T' )
cmds.menuItem( label='A' )
cmds.menuItem( label='C' )
cmds.menuItem( label='G' )
cmds.rowLayout(numberOfColumns=2)
cmds.text(label='Linker length: ')
cmds.intField('linker_length_field', value=3)
cmds.setParent('..')
cmds.button( label="Export ...", command="export_to_uniquimer_dialog()")



cmds.showWindow()

