#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sm_TM_calc_angdist.h"

int sm_TM_SNT_match(double sm_TM_RBM_matchmat[][5],
                     int L,
                     double fe_output[][2],
                     int N,
                     int sm_TM_SNT[5060][11],
                     double sm_GD_SC[5060][4],
                     double sm_TM_CP_F,
                     double pixel_size,
                     int Nth,
                     double error,
                     double sm_TM_SNT_output[][3])
{
    double fe_unmatched[N][2];
    int num_unmatched = 0;

    for (int j=0; j < N; j++)
    {
        bool already_matched = false;

        for (int i=0; i < L; i++)
        {
            if ((sm_TM_RBM_matchmat[i][2] == fe_output[j][0]) and (sm_TM_RBM_matchmat[i][3] == fe_output[j][1]))
            {
                already_matched = true;
                break;
            }
        }
        if (not already_matched)
        {
            fe_unmatched[num_unmatched][0] = fe_output[j][0];
            fe_unmatched[num_unmatched][1] = fe_output[j][1];
            num_unmatched++;
        }
    }

    int new_matched = 0, running = 0;

for (running = 0; running < num_unmatched; running++)
{
    bool done = false;
    for (int i=0; i < L; i++)
    {
        if (new_matched == Nth - L) break;
        int curr_ref_star = sm_TM_RBM_matchmat[i][4];
        if (curr_ref_star > 5060) continue;
        int j = 0;

        while (sm_TM_SNT[curr_ref_star-1][j] != 0)
        {
            /*for (int l=0; l<(Nth-L); l++)
                {cout<<"gone"<<' '<<i<<endl;
                if (sm_TM_SNT_output[l] == sm_TM_SNT[curr_ref_star-1][j]) {done = true; break;}}

            if (done) continue;*/
            double x1 = sm_TM_RBM_matchmat[i][2];
            double y1 = sm_TM_RBM_matchmat[i][3];
            double x2 = fe_unmatched[running][0];
            double y2 = fe_unmatched[running][1];

            double q1 = sm_TM_id_angdist((curr_ref_star-1), (sm_TM_SNT[curr_ref_star-1][j]-1), sm_GD_SC);
            double q2 = sm_TM_centroid_angdist(x1, y1, x2, y2, sm_TM_CP_F, pixel_size);
            double q = q1 - q2;
            //cout<<q1<<' '<<q2<<' '<<error<<endl;

            if (q<error and q>(-error))
            {
                sm_TM_SNT_output[new_matched][2] = sm_TM_SNT[curr_ref_star-1][j];
                sm_TM_SNT_output[new_matched][0] = fe_unmatched[running][0];
                sm_TM_SNT_output[new_matched][1] = fe_unmatched[running][1];

                new_matched++;
                done = true;
            }
            j++;
            if (done) break;

        }
        if (done) break;
    }
}

return new_matched;

}
