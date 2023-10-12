#include <stdio.h>
#include "FuzzyPID.h"
range rang = { 1000,-1000,800,-800,100,-100,0.1,-0.1,0.01,-0.01 };
Error error = { 0,0,0,0 };
int main()
{
    FuzzyPID myfuzzypid;
    FuzzyPID_Init(&myfuzzypid);
    float Target =1000;//目标值
    float actual = 0;//实际值
    for (int i = 0; i < 50; i++)
    {
        float u;
        u = FuzzyPIDcontroller(&myfuzzypid, &rang, &error, Target, actual);
        actual += u;
        printf("i:%d\tTarget:%f\tActual:%f\t\n",i,Target,actual);
        //std::cout << "i:" << i << "\t" << "Target:" << Target << "\t" << "Actual:" << actual << std::endl;
    }
}
