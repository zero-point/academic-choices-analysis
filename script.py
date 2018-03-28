#__future__ import print_function
import sys
import json
import math
import time
import plotly.plotly as py
import plotly.graph_objs as go
from collections import defaultdict

#time1 = time.time() # saving the start time for later

## Plotting the graph

def graph(ent,ent2,ent3,n):

        trace1 = go.Scatter(
            x=range(1,len(ent)+1),
            y=ent,
	    name='total',
	    mode='markers'
        )
        trace2 = go.Scatter(
            x=range(1,len(ent2)+1),
            y=ent2,
	    name='first event pv contribution',
	    mode='markers'
        )
        trace3 = go.Scatter(
            x=range(1,len(ent3)+1),
            y=ent3,
	    name='page view contribution',
	    mode='markers'
        )

        data = [trace1,trace3]
	layout = go.Layout(
	    title=str(n)+'gram',
	    xaxis=dict(
		title='sequences',
		titlefont=dict(
		    family='Courier New, monospace',
		    size=18,
		    color='#7f7f7f'
		)
	    ),
	    yaxis=dict(
		title='entropy',
		titlefont=dict(
		    family='Courier New, monospace',
		    size=18,
		    color='#7f7f7f'
		)
	    )
	)    
	fig = go.Figure(data=data, layout=layout)
	plot_url = py.plot(fig, filename='project'+str(n))

## opening file, reading contents

def readFile(filename):

        with open(filename) as json_file:
            json_data = json.load(json_file)
        return json_data

## calculating the probability of an event

def prob(e,a):	# c/t the no of occurences of e in a / no of total events in a

	c = 0.0
	t = 0.0
	for i in a:
		for j in i:
			t += 1
			if j==e:
				c+=1
	return c/t

def single(e,a):	# c/t the no of occurences of e in a / no of total events in a

	c = 0.0
	t = 0.0
	for i in a:
		for j in i:
			t += len(j)
			c+=j.count(e)
	return c/t

def total_events(a):
	c = 0.0
	t = 0.0
	for i in a:
		for j in i:
			t += 1
	return t

def count_event(e,a):
	c = 0.0
	for i in a:
		for j in i:
			if str(j)==str(e):
			#	print i
				c+=1
	return c


# calculate the entropy of that seq in relation to the whole
# log of geo mean = 1/N (product of log of prob of x_i)

def find_prob(a,lex,n):
	p = dict()
	for e in lex:
		p[e] = float(prob(e,a))
	return p

def entr(a,seq,lex,p):
        res = 0
        for e in lex:
		c = seq.count(e)
		if c!=0:
			res += c*math.log10(p[e])
		else:
			res += 0
	res = res*(1.0/len(seq))
        return res


def generate(access,n):

	lex = list()
	a = list()
	str_access = list()
	
	for i in access:
		s = ""
		for j in i:
			s += j
		str_access.append(s)

	for seq in str_access:
		ac = list()
		s = ""
		se = ""
		for i in range(1,n):
			se += "f" 
		se += seq
		for k in range(0,len(se)-(n-1)):
			ac.append(se[k:k+n])	 
			if se[k:k+n] not in lex:
				lex.append(se[k:k+n])
		a.append(ac)
	return lex,a

#########################
##                      #
##      Main            #
##                      #
#########################

li = list() # list of events in one access
li2 = list() # list of all accesses
lex = list()
json_data = readFile("test.json")
for user_id in json_data.keys():
	for access in json_data[user_id]:
		li = list()
		for event in access:
			e = str(event[0])
			if(e=='10'): e = 'a'
			elif(e=='11'): e = 'b'
			elif(e=='12'): e = 'c'
			elif(e=='13'): e = 'd'
			li.append(e)
			if e not in lex:
				lex.append(e)
		li2.append(li)

n = int(sys.argv[1])
a = li2
if n>1:
	lex,a=generate(a,n)
p = find_prob(a,lex,n)

ent = list()
ent2 = list()
ent3 = list()
matrix = defaultdict(dict)
m = dict()

te = list()
tee = list()
d = dict()
di = dict()
t = -10
count = 0
for seq in a:
	temp = entr(a,seq,lex,p)
	if 'f0' in seq:
		temp2 = entr(a,['f0'],lex,p)
	else:
		temp2 = temp
	ent2.append(temp2)
	if '00' in seq:
		temp3 = entr(a,['00'],lex,p)
	else:
		temp3 = temp
	ent3.append(temp3)
	
	ent.append(temp)
	m[temp] = temp3	

	if len(seq) == 1:
		count += 1

	if temp > t: 
		t = temp
		te = seq
		if len(tee) < 9:
			tee.append(te)
			d[t] = te

ent = sorted(ent, key=float, reverse=True)
ent2 = sorted(ent2, key=float, reverse=True)
ent3 = sorted(ent3, key=float, reverse=True)


graph(ent,ent2,ent3,n)























