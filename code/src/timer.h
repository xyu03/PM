
#include <ctime>
#include <cstdio>
using namespace std;

double begin_time;

void GetBeginTime()
{
    begin_time=clock()/CLOCKS_PER_SEC;
}

double GetTime()
{
    return clock()/CLOCKS_PER_SEC;
}