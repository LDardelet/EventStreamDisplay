from event import Event
import tools 

import numpy as np
import os
import pprofile
import time
import select
import socket
import cPickle

mainAddress = ("localhost", 54242)

class Player:
    def __init__(self):
        self.Socket = None
        self.Streams = {}
        self.StreamsHistory = []

    def PlayStream(self, StreamName = None, Stream = None, verboseRatio = 100000, filter_same_TS = False, resume = False):
        if not resume:
            self.nEvent = 0

        if (Stream == None and StreamName == None) or (Stream != None and StreamName != None):
	    return None

        DisplayUp = tools.IsDisplayUp()
        if not DisplayUp:
            return None
        
        if StreamName != None and StreamName not in self.Streams.keys():
            self.Streams[StreamName] = tools.load_data(StreamName, events_restriction = [0, 10**7], filter_same_TS = filter_same_TS)
            self.StreamsHistory += [StreamName]
	    Stream = self.Streams[StreamName][self.nEvent:]
	elif Stream != None:
	    StreamName = 'Imported'
	    
        if self.Socket == None:
            self.Socket = tools.GetDisplaySocket()
        else:
            tools.CleanMapForStream(self.Socket)

        tools.SendStreamData('Player',StreamName, self.Socket)
        MainUDP = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)

        starttime = time.time()
        if Stream[0].__class__.__name__ == 'Event':
            startTs = Stream[0].timestamp
        NEvents = len(Stream)
        
        for NextEvent in Stream:
            if NextEvent.__class__.__name__ == 'Event':
                self.nEvent += 1
                
                if self.nEvent%verboseRatio == 0:
                    print "New event : {1}/{2}. Polarity {0} at {3}".format(NextEvent.polarity, self.nEvent, NEvents, NextEvent.location)
                if sys.stdin in select.select([sys.stdin], [], [], 0)[0]:
                    if 'q' in sys.stdin.readline():
                        print "Closed main loop at event {0}".format(self.nEvent)
                        break
                while time.time() - starttime < NextEvent.timestamp - startTs:
                    time.sleep(0.00001)
                tools.SendEvent(NextEvent, MainUDP, mainAddress, self.Socket)
            elif NextEvent.__class__.__name__ == 'Segment':
                self.nEvent += 1
                            
                if self.nEvent%verboseRatio == 0:                
                    print "New Segment between : {0} and {1}".format(NextEvent.Points[0,:], NextEvent.Points[1,:])+NextEvent.Active*" activated" + (1-NextEvent.Active)* " deactivated"            
                if sys.stdin in select.select([sys.stdin], [], [], 0)[0]:                
                    if 'q' in sys.stdin.readline():                    
                        print "Closed main loop at event {0}".format(self.nEvent)                    
                        break            
                time.sleep(0.0001)            
                tools.SendSegment(NextEvent, MainUDP, mainAddress, self.Socket)

    def ResumeStreamPlay(self, Stream, verboseRatio = 100000):
        self.PlayStream(Stream, verboseRatio, resume = True)

