from event import Event

import numpy as np
import random
import pprofile
import sys
import pickle
import cPickle
import socket
from struct import unpack, pack

def IsDisplayUp(questionPort = 54243, responsePort = 54244, address = "localhost"):
    ResponseUDP = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
    listen_addr = ("",responsePort)
    ResponseUDP.bind(listen_addr)
    id_random = random.randint(100000,200000)
    QuestionDict = {'id': id_random, 'command':'isup'}
    QuestionUDP = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
    QuestionAddress = (address, questionPort)
    QuestionUDP.sendto(cPickle.dumps(QuestionDict),QuestionAddress)
    QuestionUDP.close()
    ResponseUDP.settimeout(1.)
    try:
        data_raw, addr = ResponseUDP.recvfrom(1064)
        data = cPickle.loads(data_raw)
        if data['id'] == id_random and data['answer']:
            return True
        else:
            return False
    except:
        print("No answer, display is down")
        return False

def KillDisplay(mainPort = 54242, questionPort = 54243, address = "localhost"):
    Question = "kill#"
    QuestionUDP = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
    QuestionAddress = (address, questionPort)
    QuestionUDP.sendto(Question,QuestionAddress)
    QuestionUDP.close()

    Main = "kill#"
    MainUDP = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
    MainAddress = (address, mainPort)
    MainUDP.sendto(Main,MainAddress)
    MainUDP.close()

def DestroySocket(Socket, questionPort = 54243, responsePort = 54244, address = "localhost"):
    if Socket == None:
        return None
    ResponseUDP = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    listen_addr = ("",responsePort)
    ResponseUDP.bind(listen_addr)

    id_random = random.randint(100000,200000)
    QuestionDict = {'id': id_random, 'socket': Socket, 'command':'destroysocket'}
    QuestionUDP = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
    QuestionAddress = (address, questionPort)
    QuestionUDP.sendto(cPickle.dumps(QuestionDict),QuestionAddress)
    QuestionUDP.close()
    ResponseUDP.settimeout(1.)

    try:
        data_raw, addr = ResponseUDP.recvfrom(1064)
        data = cPickle.loads(data_raw)
        if data['id'] == id_random and data['answer'] == 'socketdestroyed':
            print("Destroyed socket {0}".format(Socket))
        else:
            print("Could not destroy socket {0}".format(Socket))
            return None
    except:
        print("Display seems down (DestroySocket)")
    ResponseUDP.close()

def GetDisplaySocket(Socket = None, questionPort = 54243, responsePort = 54244, address = "localhost"):
    ResponseUDP = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    listen_addr = ("",responsePort)
    ResponseUDP.bind(listen_addr)

    id_random = random.randint(100000,200000)
    QuestionDict = {'id': id_random}
    if Socket == None:
        QuestionDict['command'] = "asksocket"
    else:
        QuestionDict['command'] = "askspecificsocket"
        QuestionDict['socket'] = Socket
    QuestionUDP = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
    QuestionAddress = (address, questionPort)
    QuestionUDP.sendto(cPickle.dumps(QuestionDict),QuestionAddress)
    QuestionUDP.close()
    ResponseUDP.settimeout(1.)

    try:
        data_raw, addr = ResponseUDP.recvfrom(1064)
        data = cPickle.loads(data_raw)
        if data['id'] == id_random:
            if data['answer'] == 'socketexists':
                print("Socket refused")
                return None
            else: 
                Socket = data['answer']
                print("Got socket {0}".format(Socket))
                return Socket
        else:
            print("Socket refused")
            return None
    except:
        print("Display seems down (GetDisplaySocket)")
    ResponseUDP.close()

def CleanMapForStream(Socket, questionPort = 54243, responsePort = 54244, address = "localhost"):
    ResponseUDP = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    listen_addr = ("",responsePort)
    ResponseUDP.bind(listen_addr)

    id_random = random.randint(100000,200000)
    QuestionDict = {'id': id_random, 'socket': Socket, 'command':'cleansocket'}
    QuestionUDP = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
    QuestionAddress = (address, questionPort)
    QuestionUDP.sendto(cPickle.dumps(QuestionDict),QuestionAddress)
    QuestionUDP.close()
    ResponseUDP.settimeout(1.)

    try:
        data_raw, addr = ResponseUDP.recvfrom(1064)
        data = cPickle.loads(data_raw)
        if data['id'] == id_random and data['answer'] == 'socketcleansed':
            print("Cleansed")
        else:
            print("Could not clean map")
    except:
        print("Display seems down (CleanMapForStream)")
    ResponseUDP.close()

def SendStreamData(Infos1, Infos2, Socket, questionPort = 54243, responsePort = 54244, address = "localhost"):
    ResponseUDP = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    listen_addr = ("",responsePort)
    ResponseUDP.bind(listen_addr)

    id_random = random.randint(100000,200000)
    QuestionDict = {'id': id_random, 'socket': Socket, 'infosline1': Infos1, 'infosline2': Infos2, 'command':'socketdata'}
    QuestionUDP = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)
    QuestionAddress = (address, questionPort)
    QuestionUDP.sendto(cPickle.dumps(QuestionDict),QuestionAddress)
    QuestionUDP.close()
    ResponseUDP.settimeout(1.)

    try:
        data_raw, addr = ResponseUDP.recvfrom(1064)
        data = cPickle.loads(data_raw)
        if data['id'] == id_random and data['answer'] == 'datareceived':
            print("Data transmitted")
        else:
            print("Could not transmit data")
    except:
        print("Display seems down (SendStreamData)")
    ResponseUDP.close()


def SendEvent(ev, Connexion, Address, Socket):
    ev.socket = Socket
    data = pickle.dumps(ev)
    Connexion.sendto(data, Address)

def SendSegment(segment, Connexion, Address, Socket):
    segment.socket = Socket
    data = pickle.dumps(segment)
    Connexion.sendto(data, Address)
    
def StatFunction(function, *args):
    prof = pprofile.Profile()

    with prof():
        try:
            function(*args)
        except KeyboardInterrupt:
            print("Function ended, writing stats in file")
    
    orig_stdout = sys.stdout
    f = file('stats.txt', 'w')
    sys.stdout = f
    prof.print_stats()
    sys.stdout = orig_stdout
    f.close()

def load_data(filename, y_invert = True, listeningPolas = [0,1], header_size = 0, geometry_restriction = [[10,230], [10, 230]], time_restriction = [0., np.inf], events_restriction = [0, 5000000]):
    if '.dat' in filename:
        all_ts, coords, all_p, rem = readATIS_td(filename)
        all_x = coords[:,0]
        all_y = coords[:,1]
        y_max = geometry_restriction[1][1]
        del coords
    else:
        f = open(filename, 'rb')
        raw_data = np.fromfile(f, dtype=np.uint8)[header_size:]
        f.close()
        raw_data = np.uint32(raw_data)
        all_y = raw_data[1::5]
        y_max = max(all_y)
        all_x = raw_data[0::5]
        all_p = (raw_data[2::5] & 128) >> 7
        all_ts = ((raw_data[2::5] & 127) << 16) | (raw_data[3::5] << 8) | (raw_data[4::5])

    print("Read file. Seems to contain {0} events until ts = {1}s".format(all_x.shape[0], all_ts[-2]*10**-6))
    print("")
    print("Geometry :")
    print("x goes from {0} to {1} (restricted to {2}-{3})".format(all_x.min(), all_x.max(), geometry_restriction[0][0], geometry_restriction[0][1]))
    print("y goes from {0} to {1} (restricted to {2}-{3})".format(all_y.min(), all_y.max(), geometry_restriction[1][0], geometry_restriction[1][1]))
    Events = []
    for nEvent in range(all_x.shape[0]):
    #for x,y,p,ts in zip(all_x, all_y, all_p, all_ts):
        x,y,p,ts = all_x[nEvent], all_y[nEvent], all_p[nEvent], all_ts[nEvent]
        if p in listeningPolas:
            if y_invert:
                y_final = y_max - y
            else:
                y_final = y
            t = ts*(10**-6)
            if events_restriction[0] <= nEvent <= events_restriction[1] and time_restriction[0] <= t <= time_restriction[1]:
                X = int(x)
                if geometry_restriction[0][0] <= X <= geometry_restriction[0][1]:
                    Y = int(y_final)
                    if geometry_restriction[1][0] <= Y <= geometry_restriction[1][1]:
                        Events += [Event(float(t), np.array([X-geometry_restriction[0][0],Y-geometry_restriction[1][0]]), int(listeningPolas.index(p)))]
                        if nEvent%10000 == 0:
                            print(Events[-1].timestamp, nEvent)
            elif nEvent> events_restriction[1] or t > time_restriction[1]:
                break
    return Events

def load_data_csv(filename, listeningPolas = [0,1], eventType = [0], geometry_restriction = [[10,230], [10, 230]], time_restriction = [1., np.inf]):
    Events = []

    with open(filename) as csvfile:
        for raw_line in csvfile:
            line = raw_line.split(',')
            if int(line[3].strip('\t').strip('\n')) in eventType:
                P = int(line[4].strip('\t').strip('\n'))
                if P in listeningPolas:
                    t = float(line[2].strip('\t').strip('\n'))*(10**-6)
                    if time_restriction[0] <= t <= time_restriction[1]:
                        X = int(line[0].strip('\t').strip('\n'))
                        if geometry_restriction[0][0] <= X <= geometry_restriction[0][1]:
                            Y = int(line[1].strip('\t').strip('\n'))
                            if geometry_restriction[1][0] <= Y <= geometry_restriction[1][1]:
                                Events += [Event(t, np.array([X-geometry_restriction[0][0],Y-geometry_restriction[1][0]]), P)]
                                if random.random() < 0.001:
                                    print Events[-1].timestamp
    return Events

def load_StreamClass(filename):
    try:
        return filename.split('/')[-2]
    except:
        return 'unknown'

def CreateTimeSurface(History, ultimate_time, tau, time_forget, min_events_number_in_TS):
    if (History > ultimate_time-time_forget).sum() >= min_events_number_in_TS:
        TimeSurface = np.e**(-(ultimate_time-History)/tau)
        return TimeSurface
    else:
        return None

def Insert(Scene, Patch, x, y, getEventsDiff = False):
	if getEventsDiff:
		InitialScene = np.array(Scene)
	RScene = (Scene.shape[0]-1)/2
	RPatch = (Patch.shape[0]-1)/2
	minpatchx = max(0,RPatch-x)
	maxpatchx = 2*RPatch+1 - max(0,x - 2*RScene + RPatch)
	minpatchy = max(0,RPatch-y)
	maxpatchy = 2*RPatch+1 - max(0,y - 2*RScene + RPatch)
	ScenePart = np.array(Scene[max(0,x-RPatch):x+RPatch+1,max(0,y-RPatch):y+RPatch+1,:])
	Scene[max(0,x-RPatch):x+RPatch+1,max(0,y-RPatch):y+RPatch+1,:] = np.maximum(Patch[minpatchx:maxpatchx, minpatchy:maxpatchy, :], ScenePart)
	if getEventsDiff:
		return Scene, np.where(Scene != InitialScene)
	else:
		return Scene
	
def ComputeFlowMove(F, Neighbours):
    F = np.array(F)/np.linalg.norm(F)
    best_dis = 0
    best_index = 0
    for nN in range(len(Neighbours)):
        d = (F*Neighbours[nN]).sum()
        if d > best_dis:
            best_index = nN
            best_dis = d
    return best_index, best_dis

def GetOrthogonalNormalizedVector(Vector):
    Vector /= np.linalg.norm(Vector)
    return np.array([-Vector[1], Vector[0]])

def readATIS_td(file_name, orig_at_zero = True, drop_negative_dt = True, verbose =
True):

    # This one read _td.dat files generated by kAER
    if verbose:
        print('Reading _td dat file... (' + file_name + ')')
    file = open(file_name,'rb')

    header = False
    while peek(file) == b'%':
        file.readline()
        header = True
    if header:
        ev_type = unpack('B',file.read(1))[0]
        ev_size = unpack('B',file.read(1))[0]
        if verbose:
            print('> Header exists. Event type is ' + str(ev_type) + ', event size is ' + str(ev_size))
        if ev_size != 8:
            print('Wrong event size. Aborting.')
            return -1, -1, -1, -1
    else: # set default ev type and size
        if verbose:
            print('> No header. Setting default event type and size.')
        ev_size = 8
        ev_type = 0

    # Compute number of events in the file
    start = file.tell()
    file.seek(0,2)
    stop = file.tell()
    file.seek(start)

    Nevents = int( (stop-start)/ev_size )
    dNEvents = Nevents/100
    if verbose:
        print("> The file contains %d events." %Nevents)

    # store read data
    timestamps = np.zeros(Nevents, dtype = int)
    polarities = np.zeros(Nevents, dtype = int)
    coords = np.zeros((Nevents, 2), dtype = int)

    for i in np.arange(0, int(Nevents)):

        event = unpack('Q',file.read(8))
        ts = event[0] & 0x00000000FFFFFFFF
        # padding = event[0] & 0xFFFC000000000000
        pol = (event[0] & 0x0002000000000000) >> 49
        y = (event[0] & 0x0001FE0000000000) >> 41
        x = (event[0] & 0x000001FF00000000) >> 32

        timestamps[i] = ts
        polarities[i] = pol
        coords[i, 0] = x
        coords[i, 1] = y

        if verbose and i%dNEvents == 0:
            sys.stdout.write("> "+str(i/dNEvents)+"% \r")
            sys.stdout.flush()

    file.close()

    if verbose:
        print("> Sequence duration: %ds" %float((timestamps[-1] - timestamps[0]) / 1e6))

    #check for negative timestamps
    for ts in timestamps:
        if ts < 0:
            print('prout negative delta-ts')

    if orig_at_zero:
        timestamps = timestamps - timestamps[0]

    drop_sum = 0
    if drop_negative_dt:
        if verbose:
            print('> Looking for negative dts...')
        # first check if negative TS differences
        just_dropped = True
        nPasses = 0
        while just_dropped:
            nPasses += 1
            index_neg = []
            just_dropped = False
            ii = 0
            while ii < (timestamps.size - 1):
                dt = timestamps[ii+1] - timestamps[ii]
                if dt <= 0:  # alors ts en ii+1 plus petit que ii
                    index_neg += [ii+1]
                    ii += 1
                    just_dropped = True
                if verbose and ii%dNEvents == 0:
                    sys.stdout.write("> "+str(ii/dNEvents)+"% (pass "+str(nPasses)+") \r")
                    sys.stdout.flush()
                ii += 1
            if len(index_neg) > 0:
                drop_sum += len(index_neg)
                index_neg = np.array(index_neg)
                timestamps = np.delete(timestamps, index_neg)
                polarities = np.delete(polarities, index_neg)
                coords = np.delete(coords, index_neg, axis = 0)
                if verbose:
                    print('> Removed {0} events in {1} passes.'.format(drop_sum, nPasses))
        removed_events = drop_sum
    else:
        removed_events = -1

    return timestamps, coords, polarities, removed_events

def peek(f, length=1):
    pos = f.tell()
    data = f.read(length)
    f.seek(pos)
    return data

def CompareStreams(S1, S2, geometry, deltaTMax):
	TS1 = 0
	TS2 = 0
	
	Sum1 = 0
	Sum2 = 0
	
	TSList1 = [[[[] for p in range(geometry[2])] for y in range(geometry[1])] for x in range(geometry[0])]
	TSList2 = [[[[] for p in range(geometry[2])] for y in range(geometry[1])] for x in range(geometry[0])]
	
	for ev1 in S1:
		TSList1[ev1.location[0]][ev1.location[1]][ev1.polarity] += [ev1.timestamp]
	for ev2 in S2:
		TSList2[ev2.location[0]][ev2.location[1]][ev2.polarity] += [ev2.timestamp]
	
	for x in range(geometry[0]):
		sys.stdout.write("Treating line {0}/{1} \r".format(x, geometry[0]-1))
		sys.stdout.flush()
		for y in range(geometry[1]):
			for p in range(geometry[2]):
				Pixel1 = list(TSList1[x][y][p])
				Pixel2 = list(TSList2[x][y][p])
				L1 = len(Pixel1)
				L2 = len(Pixel2)
				if L1 != 0 and L2 != 0:
					nts = 0
					for ts in Pixel1:
						nts += 1
						diff2 = abs(np.array(Pixel2) - ts)
						index2 = diff2.argmin()
						if diff2[index2] > deltaTMax:
							Sum1 += deltaTMax
						else:
							Sum1 += diff2[index2]
							Pixel2.pop(index2)
						if len(Pixel2) == 0:
							Sum1 += (L1-nts)*deltaTMax
							break
				
					Pixel1 = list(TSList1[x][y][p])
					Pixel2 = list(TSList2[x][y][p])
				
					nts = 0
					for ts in Pixel2:
						nts += 1
						diff1 = abs(np.array(Pixel1) - ts)
						index1 = diff1.argmin()
						if diff1[index1] > deltaTMax:
							Sum2 += deltaTMax
						else:
							Sum2 += diff1[index1]
							Pixel1.pop(index1)
						if len(Pixel1) == 0:
							Sum2 += (L2-nts)*deltaTMax
							break
				else:
					Sum1 += L1*deltaTMax
					Sum2 += L2*deltaTMax
	print ""
	print "Sum1 : {0}".format(Sum1)
	print "Sum2 : {0}".format(Sum2)
	return (Sum1+Sum2)/2
