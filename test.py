from math import cos, pi
tot_day = 8.5

hour=-6
res=0.5
tot=0.0;
while hour < 6:
    rtime=hour*pi/12.0;
    y=tot_day*cos(rtime)*(res*pi/12.0)/2;
    tot+=y;
    hour+=res;

print tot, tot_day
