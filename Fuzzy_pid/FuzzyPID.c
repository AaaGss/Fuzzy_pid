#include "FuzzyPID.h"
#define NB -3
#define NM - 2
#define NS - 1
#define ZO 0
#define PS 1
#define PM 2
#define PB 3
int  Kp_rule[7][7] = { {PB,PB,PM,PM,PS,ZO,ZO},     //kp规则表
                               {PB,PB,PM,PS,PS,ZO,NS},
                               {PM,PM,PM,PS,ZO,NS,NS},
                               {PM,PM,PS,ZO,NS,NM,NM},
                               {PS,PS,ZO,NS,NS,NM,NM},
                               {PS,ZO,NS,NM,NM,NM,NB},
                               {ZO,ZO,NM,NM,NM,NB,NB} };

int  Ki_rule[7][7] = { {NB,NB,NM,NM,NS,ZO,ZO},     //ki规则表
                            {NB,NB,NM,NS,NS,ZO,ZO},
                            {NB,NM,NS,NS,ZO,PS,PS},
                            {NM,NM,NS,ZO,PS,PM,PM},
                            {NM,NS,ZO,PS,PS,PM,PB},
                            {ZO,ZO,PS,PS,PM,PB,PB},
                            {ZO,ZO,PS,PM,PM,PB,PB} };

int  Kd_rule[7][7] = { {PS,NS,NB,NB,NB,NM,PS},    //kd规则表
                            {PS,NS,NB,NM,NM,NS,ZO},
                            {ZO,NS,NM,NM,NS,NS,ZO},
                            {ZO,NS,NS,NS,NS,NS,ZO},
                            {ZO,ZO,ZO,ZO,ZO,ZO,ZO},
                            {PB,NS,PS,PS,PS,PS,PB},
                            {PB,PM,PM,PM,PS,PS,PB} };

int  Fuzzy_rule[7][7] = { {PB,PB,PB,PB,PM,ZO,ZO},
                               {PB,PB,PB,PM,PM,ZO,ZO},
                               {PB,PM,PM,PS,ZO,NS,NM},
                               {PM,PM,PS,ZO,NS,NM,NM},
                               {PS,PS,ZO,NM,NM,NM,NB},
                               {ZO,ZO,ZO,NM,NB,NB,NB},
                               {ZO,NS,NB,NB,NB,NB,NB} };
float values[7] = { -3,-2,-1,0,1,2,3 }; //输入e的隶属值
void FuzzyPID_Init(FuzzyPID* pid)  //构造函数
{
    int i, j;
    pid->num_area = 8;
    pid->kp = 0;
    pid->ki = 0;
    pid->kd = 0;
    pid->fuzzy_output = 0;
    pid->qdetail_kp = 0;
    pid->qdetail_ki = 0;
    pid->qdetail_kd = 0;
    pid->qfuzzy_output = 0;
    pid->errosum = 0;
    for ( i = 0; i < 7; i++)
    {
        for ( j = 0; j < 7; j++)
        {
            pid->Kp_rule_list[i][j] = Kp_rule[i][j];
            pid->Ki_rule_list[i][j] = Ki_rule[i][j];
            pid->Kd_rule_list[i][j] = Kd_rule[i][j];
            pid->Fuzzy_rule_list[i][j] = Fuzzy_rule[i][j];
        }
    }
    for ( i = 0; i < 7; i++)
    {
        pid->e_membership_values[i] = values[i];
        pid->ec_membership_values[i] = values[i];
        pid->kp_menbership_values[i] = values[i];
        pid->ki_menbership_values[i] = values[i];
        pid->kd_menbership_values[i] = values[i];
        pid->fuzzyoutput_menbership_values[i] = values[i];
        pid->gradSums[i] = 0;
        pid->KpgradSums[i] = 0;
        pid->KigradSums[i] = 0;
        pid->KdgradSums[i] = 0;
    }
}




//输入e与de/dt隶属度计算函数///
void Get_grad_membership(FuzzyPID* pid,float erro, float erro_c)
{
    int i;
    //当误差在这个范围
    if (erro > pid->e_membership_values[0] && erro < pid->e_membership_values[6])
    {

        //6个区域
        for ( i = 0; i < pid->num_area - 2; i++)
        {
            //如果误差在区间区域内
            if (erro >= pid->e_membership_values[i] && erro <= pid->e_membership_values[i + 1])
            {
                //e的隶属度
                //PM
                pid->e_gradmembership[0] = -(erro - pid->e_membership_values[i + 1]) / (pid->e_membership_values[i + 1] - pid->e_membership_values[i]);
                //PB
                pid->e_gradmembership[1] = 1 + (erro - pid->e_membership_values[i + 1]) / (pid->e_membership_values[i + 1] - pid->e_membership_values[i]);
                //记录是在哪两个区间内
                pid->e_grad_index[0] = i;
                pid->e_grad_index[1] = i + 1;
                break;
            }
        }
    }
    else
    {
        //如果误差的止小于等于论域的最小值
        if (erro <= pid->e_membership_values[0])
        {
            pid->e_gradmembership[0] = 1;
            pid->e_gradmembership[1] = 0;
            pid->e_grad_index[0] = 0;
            pid->e_grad_index[1] = -1;
        }//超出范围了
        else if (erro >= pid->e_membership_values[6])
        {
            pid->e_gradmembership[0] = 1;
            pid->e_gradmembership[1] = 0;
            pid->e_grad_index[0] = 6;
            pid->e_grad_index[1] = -1;
        }
    }
    //误差的微分
    if (erro_c > pid->ec_membership_values[0] && erro_c < pid->ec_membership_values[6])
    {
        for ( i = 0; i < pid->num_area - 2; i++)
        {
            if (erro_c >= pid->ec_membership_values[i] && erro_c <= pid->ec_membership_values[i + 1])
            {
                pid->ec_gradmembership[0] = -(erro_c - pid->ec_membership_values[i + 1]) / (pid->ec_membership_values[i + 1] - pid->ec_membership_values[i]);
                pid->ec_gradmembership[1] = 1 + (erro_c - pid->ec_membership_values[i + 1]) / (pid->ec_membership_values[i + 1] - pid->ec_membership_values[i]);
                pid->ec_grad_index[0] = i;
                pid->ec_grad_index[1] = i + 1;
                break;
            }
        }
    }
    else
    {
        if (erro_c <= pid->ec_membership_values[0])
        {
            pid->ec_gradmembership[0] = 1;
            pid->ec_gradmembership[1] = 0;
            pid->ec_grad_index[0] = 0;
            pid->ec_grad_index[1] = -1;
        }
        else if (erro_c >= pid->ec_membership_values[6])
        {
            pid->ec_gradmembership[0] = 1;
            pid->ec_gradmembership[1] = 0;
            pid->ec_grad_index[0] = 6;
            pid->ec_grad_index[1] = -1;
        }
    }

}

// //获取输出增量kp, ki, kd的总隶属度 /
void GetSumGrad(FuzzyPID* pid)
{
    int i,j;
    //划分八个区域
    for ( i = 0; i <= pid->num_area - 1; i++)
    {
        pid->KpgradSums[i] = 0;
        pid->KigradSums[i] = 0;
        pid->KdgradSums[i] = 0;
        //把PID的各个隶属值清零
    }
    for ( i = 0; i < 2; i++)//循环两次
    {
        if (pid->e_grad_index[i] == -1)//误差有没有爆表
        {
            continue;
        }
        for ( j = 0; j < 2; j++)//
        {
            if (pid->ec_grad_index[j] != -1)//误差的微分有没有爆表
            {
                int indexKp = pid->Kp_rule_list[pid->e_grad_index[i]][pid->ec_grad_index[j]] + 3;
                int indexKi = pid->Ki_rule_list[pid->e_grad_index[i]][pid->ec_grad_index[j]] + 3;
                int indexKd = pid->Kd_rule_list[pid->e_grad_index[i]][pid->ec_grad_index[j]] + 3;
                //gradSums[index] = gradSums[index] + (e_gradmembership[i] * ec_gradmembership[j])* Kp_rule_list[e_grad_index[i]][ec_grad_index[j]];
                pid->KpgradSums[indexKp] = pid->KpgradSums[indexKp] + (pid->e_gradmembership[i] * pid->ec_gradmembership[j]);
                pid->KigradSums[indexKi] = pid->KigradSums[indexKi] + (pid->e_gradmembership[i] * pid->ec_gradmembership[j]);
                pid->KdgradSums[indexKd] = pid->KdgradSums[indexKd] + (pid->e_gradmembership[i] * pid->ec_gradmembership[j]);
            }
            else
            {
                continue;
            }

        }
    }

}

//计算输出增量kp, kd, ki对应论域值//
void GetOUT(FuzzyPID* pid)
{
    int i;
    for ( i = 0; i < pid->num_area - 1; i++)
    {
        pid->qdetail_kp +=pid->kp_menbership_values[i] * pid->KpgradSums[i];
        pid->qdetail_ki += pid->ki_menbership_values[i] * pid->KigradSums[i];
        pid->qdetail_kd += pid->kd_menbership_values[i] * pid->KdgradSums[i];
    }
}

//模糊PID控制实现函数/
float FuzzyPIDcontroller(FuzzyPID* pid, range* rang, Error* error, float Target, float actual)
{
    
    error->erro_ppre = error->erro_pre;
    error->erro_pre = error->erro;
    error->erro = Target - actual;
    error->erro_c = error->erro - error->erro_pre;
    pid->errosum += error->erro;
    //Arear_dipart(e_max, e_min, ec_max, ec_min, kp_max, kp_min,ki_max,ki_min,kd_max,kd_min);
    pid->qerro = Quantization(rang->e_max, rang->e_min, error->erro);//区间映射
    pid->qerro_c = Quantization(rang->ec_max, rang->ec_min, error->erro_c);//区间映射
    //把他们缩小到0123范围内
    Get_grad_membership(pid,pid->qerro, pid->qerro_c);
    //获取输出增量kp, ki, kd的总隶属度
    GetSumGrad(pid);
    //计算输出增量kp, kd, ki对应论域值//
    GetOUT(pid);
    pid->detail_kp = Inverse_quantization(rang->kp_max, rang->kp_min, pid->qdetail_kp);
    pid->detail_ki = Inverse_quantization(rang->ki_max, rang->ki_min, pid->qdetail_ki);
    pid->detail_kd = Inverse_quantization(rang->kd_max, rang->kd_min, pid->qdetail_kd);
    pid->qdetail_kd = 0;
    pid->qdetail_ki = 0;
    pid->qdetail_kp = 0;
    /*if (qdetail_kp >= kp_max)
        qdetail_kp = kp_max;
    else if (qdetail_kp <= kp_min)
        qdetail_kp = kp_min;
    if (qdetail_ki >= ki_max)
        qdetail_ki = ki_max;
    else if (qdetail_ki <= ki_min)
        qdetail_ki = ki_min;
    if (qdetail_kd >= kd_max)
        qdetail_kd = kd_max;
    else if (qdetail_kd <= kd_min)
        qdetail_kd = kd_min;*/
    pid->kp = pid->kp + pid->detail_kp;
    pid->ki = pid->ki + pid->detail_ki;
    pid->kd  =pid->kd + pid->detail_kd;
    //确定范围
    if (pid->kp < 0)
        pid->kp = 0;
    if (pid->ki < 0)
        pid->ki = 0;
    if (pid->kd < 0)
        pid->kd = 0;
    pid->detail_kp = 0;
    pid->detail_ki = 0;
    pid->detail_kd = 0;
    //增量式PID
    float output = pid->kp * (error->erro - error->erro_pre) + pid->ki * error->erro + pid->kd * (error->erro - 2 * error->erro_pre + error->erro_ppre);
    return output;
}

///区间映射函数///
float Quantization(float maximum, float minimum, float x)
{
    float qvalues = 6.0 * (x - minimum) / (maximum - minimum) - 3;
    //float qvalues=6.0*()
    return qvalues;

    //qvalues[1] = 3.0 * ecerro / (maximum - minimum);
}

//反区间映射函数
float Inverse_quantization(float maximum, float minimum, float qvalues)
{
    float x = (maximum - minimum) * (qvalues + 3) / 6 + minimum;
    return x;
}
