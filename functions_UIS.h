#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "constants_UIS.h"
#include "sm_K_vec_arr.h"
#include "sm_GC.h"
#include "sm_TM_SNT.h"


double absoluteValue(double x)
{
    if (x < 0) return -x;
    return x;
}

//LISM functions

void bubbleSort(double arr[][3], int n)
{
    int i, j;
    long double temp[3];
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n - i - 1; j++)
        {
            if (arr[j][1] * arr[j][1] + arr[j][2] * arr[j][2] > arr[j + 1][1] * arr[j + 1][1] + arr[j + 1][2] * arr[j + 1][2])
            {
                // swap the elements
                temp[0] = arr[j][0];
                temp[1] = arr[j][1];
                temp[2] = arr[j][2];
                arr[j][0] = arr[j + 1][0];
                arr[j][1] = arr[j + 1][1];
                arr[j][2] = arr[j + 1][2];
                arr[j + 1][0] = temp[0];
                arr[j + 1][1] = temp[1];
                arr[j + 1][2] = temp[2];
            }
        }
    }
}

void sm_4_star_circulate(double sm_3D_vecs[][4], int *N_circ, int N_i)
{   
    int v;
    int j;
    int k; 

    int last = 0;
    for (j = N_i - 1; j >= 0; j--)
    {
        if (sm_3D_vecs[j][0] != -1)
        {
            last = j;
            break;
        }
    }
    double curr[4] = {sm_3D_vecs[last][0], sm_3D_vecs[last][1], sm_3D_vecs[last][2], sm_3D_vecs[last][3]};
    double var;
    for (k = 0; k < N_i; k++)
    {
        if (sm_3D_vecs[k][0] != -1)
        {
            for (v = 0; v < 4; v++)
            {
                var = sm_3D_vecs[k][v];
                sm_3D_vecs[k][v] = curr[v];
                curr[v] = var;
            }
        }
    }
    *N_circ++;
}

int already_matched(int sm_IS[][2], int indx)
{   
    int i;
    for (i = 0; i < N_GC; i++)
    {
        if (sm_IS[i][0] == indx)
        {
            return 1;
        }
    }
    return 0;
}

void sm_4_star(double four_stars[][4], double sm_3D_vecs[][4], int sm_IS[][2], double body_vecs_IS[][4], int sm_K_vec_arr[][3], int *N_match, int N_i, double q, double m, int N_is)
{   
    int i = 0;
    int j = 0;
    int k = 0;

    int SIM[N_GC][6];
    memset(SIM, 0, N_GC * sizeof(SIM[0]));

    int indx_arr[N_GC];
    memset(indx_arr,0, N_GC*sizeof(indx_arr[0]));

    int SIM_flags[N_GC];
    memset(SIM_flags, 0, N_GC*sizeof(SIM_flags[0]));


    int top_indx = 0;

    
    double p[6];
    int ct = 0;
    for (i = 0; i < 4; i++)
    {
        for (j = i + 1; j < 4; j++)
        {
            double norm1 = sqrt(four_stars[i][1] * four_stars[i][1] + four_stars[i][2] * four_stars[i][2] + four_stars[i][3] * four_stars[i][3]);
            double norm2 = sqrt(four_stars[j][1] * four_stars[j][1] + four_stars[j][2] * four_stars[j][2] + four_stars[j][3] * four_stars[j][3]);
            p[ct] = fabs(four_stars[i][1] * four_stars[j][1] + four_stars[i][2] * four_stars[j][2] + four_stars[i][3] * four_stars[j][3]) / (norm1 * norm2);
            // printf("%f ",p[ct]);   //angular distances
            ct++;
            
        }
    }
    // printf("\n");

    int checks[4][6] = {{1, 1, 1, 0, 0, 0},
                        {1, 0, 0, 1, 1, 0},
                        {0, 1, 0, 1, 0, 1},
                        {0, 0, 1, 0, 1, 1}};

    for (j = 0; j < 6; j++)
    {
        double sin_j = sqrt(1 - (p[j] * p[j]));
        int k_top = ceil((cos(DELTA) * p[j] + sin_j * sin(DELTA) - q) / m);
        int k_bot = floor((cos(DELTA) * p[j] - sin_j * sin(DELTA) - q) / m);
        // printf(" k_bot = %d, k_top = %d\n", k_bot, k_top);
        if (k_top <= 0 || k_bot >= 224792)
        {
            printf("bad values : k_bot = %d, k_top = %d\n", k_bot, k_top);
            continue;
        }
        else
        {
            if (k_top > 224792)
            {
                k_top = 224792;
            }
            if (k_bot < 0)
            {
                k_bot = 1;
            }
            int k_start = sm_K_vec_arr[k_bot - 1][2] + 1;
            int k_end = sm_K_vec_arr[k_top - 1][2];
            // printf("k_start : %d, k_end : %d\n", k_start, k_end);
            if (k_start == k_end)
            {
                SIM[sm_K_vec_arr[k_end - 1][0] - 1][j] = 1;
                indx_arr[top_indx] = sm_K_vec_arr[k_end - 1][0] - 1 ;
                top_indx ++;
            }
            else
            {
                for (i = k_start; i <= k_end; i++)
                {
                    SIM[sm_K_vec_arr[i - 1][0] - 1][j] = 1; 
                    SIM[sm_K_vec_arr[i - 1][1] - 1][j] = 1;

                    if(SIM_flags[sm_K_vec_arr[i - 1][0] - 1] == 0) {
                        SIM_flags[sm_K_vec_arr[i - 1][0] - 1] = 1;
                        indx_arr[top_indx] = sm_K_vec_arr[i - 1][0] - 1;
                        top_indx++;
                    }

                    if(SIM_flags[sm_K_vec_arr[i - 1][1] - 1] == 0) {
                        SIM_flags[sm_K_vec_arr[i - 1][1] - 1] = 1;
                        indx_arr[top_indx] = sm_K_vec_arr[i-1][1] - 1;
                        top_indx++;
                    }
                    
                }
            }
        }
    }    
  //printing SIM 
    // printf("printing SIM\n");
    // int ip = 0;
    // for(ip = 0; ip < N_GC; ip++)
    // {
    //     printf(" %d %d %d %d %d %d %d\n", ip, SIM[ip][0], SIM[ip][1], SIM[ip][2], SIM[ip][3], SIM[ip][4], SIM[ip][5]);
    // }

    for (j = 0; j < 4; j++)
    {
        int matched_rows = 0;
        int temp = 0;  

        for(i = 0; i < top_indx; i++)
        {
            SIM_flags[indx_arr[i]] = 1;
        }

        // for(i=0; i < 8876; i++)
        // {
        //     if(SIM_flags[i] == 1)
        //         printf("%d ", i);
        // }

        // printf("--------------------\n");
        // printf("\n matched stars :");
        for (i = 0; i < top_indx; i++)
        {
            k = indx_arr[i];
            if(SIM_flags[indx_arr[i]] == 1)
            {
                if (SIM[k][0] == checks[j][0] 
                && SIM[k][1] == checks[j][1] 
                && SIM[k][2] == checks[j][2] 
                && SIM[k][3] == checks[j][3] 
                && SIM[k][4] == checks[j][4] 
                && SIM[k][5] == checks[j][5])
                {
                    matched_rows++;
                    temp = k;
                    printf("%d ",temp);
                }
                SIM_flags[k] = 0;
            }
        }
        
        if (matched_rows == 1)
        {
            // printf("matched_once\n");
            int flag = already_matched(sm_IS, (int)four_stars[j][0]);
            if (flag == 0)
            {
                (*N_match)++;
                for (k = 0; k < N_i; k++)
                {
                    if ((int)sm_3D_vecs[k][0] == (int)four_stars[j][0])
                    {
                        sm_3D_vecs[k][0] = -1;
                        break;
                    }
                }
                k = 0;
                for (k = 0; k < N_GC; k++)
                {
                    if (sm_IS[k][0] == -1)
                    {
                        sm_IS[k][0] = (int)four_stars[j][0];
                        body_vecs_IS[k][0] = four_stars[j][0];
                        sm_IS[k][1] = temp + 1;
                        //printing matched starIDs
                        // printf("%d \n", sm_IS[k][1]);
                        i = 1;
                        for (i = 1; i < 4; i++){
                            body_vecs_IS[k][i] = four_stars[j][i];
                        }
                        break;
                    }
                }
            }
        }
    }    
}

void sm_gnrt_3D_vec(double sm_3D_vecs[][4], double sm_sorted_UIS[][3], int N_i)
{   int i;
    for (i = 0; i < N_i; i++)
    {
        double local = sqrt((sm_sorted_UIS[i][1] / FOCAL_LENGTH) * (sm_sorted_UIS[i][1] / FOCAL_LENGTH) + (sm_sorted_UIS[i][2] / FOCAL_LENGTH) * (sm_sorted_UIS[i][2] / FOCAL_LENGTH) + 1);
        sm_3D_vecs[i][0] = sm_sorted_UIS[i][0];
        sm_3D_vecs[i][1] = (sm_sorted_UIS[i][1] / (FOCAL_LENGTH * local));
        sm_3D_vecs[i][2] = (sm_sorted_UIS[i][2] / (FOCAL_LENGTH * local));
        sm_3D_vecs[i][3] = 1 / local;

    }
}

void sm_validate(double sm_3D_vecs[][4], int sm_IS[][2], double body_vecs_IS[][4], double sm_GC[][4], int *N_is, int N_i, double tol, double p_1, double p_2)
{
    int i, j;
    int N_new = *N_is;
    int votes[N_i];
    memset(votes, 0, N_i * sizeof(votes[0]));
    

    for (i = 0; i < N_i; i++){
        if (sm_IS[i][0] == -1)
            continue;
        for (j = i+1; j< N_i-1; j++){
            if (sm_IS[j][0] != -1){
                double d_ij = fabs(body_vecs_IS[i][1]*body_vecs_IS[j][1] + body_vecs_IS[i][2]*body_vecs_IS[j][2] + body_vecs_IS[i][3]*body_vecs_IS[j][3]);
                double d_ij_gc = fabs(sm_GC[sm_IS[i][1] - 1][1]*sm_GC[sm_IS[j][1] - 1][1] + sm_GC[sm_IS[i][1] - 1][2]*sm_GC[sm_IS[j][1] - 1][2] + sm_GC[sm_IS[i][1] - 1][3]*sm_GC[sm_IS[j][1] - 1][3]);
                if (fabs(d_ij/d_ij_gc - 1) < tol/100){
                    votes[i]++;
                    votes[j]++;
                }
                if((votes[i] < 0) || (votes[j] < 0))
                    printf("Idx i: %d, j: %d\n", i, j);
            }
        }
    }
    int N_LB = p_1*(N_new)/100;
    for (i = 0; i < N_i; i++){
        if (sm_IS[i][0] == -1)
            continue;
        if (votes[i] < N_LB){
            sm_IS[i][0] = -1;
            body_vecs_IS[i][0] = -1;
            N_new--;
        }
    }
    int N_UB = p_2*(N_new)/100;
    for (i = 0; i < N_i; i++){
        if (sm_IS[i][0] == -1)
            continue;
        if (votes[i] < N_UB){
            sm_IS[i][0] = -1;
            body_vecs_IS[i][0] = -1;
            N_new--;
        }
    }
    *N_is = N_new;
}

void LISM(double centroids_st[MAX_STARS][3], int tot_stars, double data[3][MAX_STARS], int input_ids[50], int star_ids[50], int* matched_stars){
    // Maximum number of iterations
    int N_max = tot_stars - 1;

    // N_i is specially used to pass onto other functions(remains unchanged throughout the code)
    int N_i = tot_stars;

    // N_uis is a SM variable for the total number of stars detected from the FE block(changes value in main)
    int N_uis = tot_stars;

    // Required maximum number of matched stars
    int N_th = 10;

    // Number of times control goes to sm_4_star_circulate
    int N_circ = 0;

    // Number of identified stars
    int N_is = 0;

    double tol = TOL;
    double p_1 = P1;
    double p_2 = P2;

    // Constants for using the k vector table (to be used in the 4 star matching)
    double m = (Y_MAX - Y_MIN + 2 * EPSILON) / (N_KVEC_PAIRS - 1);
    double q = Y_MIN - m - EPSILON;

    // Array for storing the Identified Stars
    int sm_IS[N_GC][2];
    // Array for storing body-frame vectors of Identified stars
    double body_vecs_IS[N_GC][4];

    //Sorting the UIS according to euclidean distance
    bubbleSort(centroids_st, N_i);

    // Initialize a block of memory to -1
    memset(sm_IS, -1, N_GC * sizeof(sm_IS[0]));
    memset(body_vecs_IS, -1, N_GC * sizeof(body_vecs_IS[0]));

    // Array to store body frame vectors of extracted stars from FE block
    double sm_3D_vecs[N_i][4];

    // Generate body frame vectors for extracted stars
    sm_gnrt_3D_vec(sm_3D_vecs, centroids_st, N_i);

    // for(int i = 0; i < N_i; i++)
    // {
    //     printf("%f %f %f %f\n", sm_3D_vecs[i][0],  sm_3D_vecs[i][1],  sm_3D_vecs[i][2],  sm_3D_vecs[i][3]);
    // }

    int i = 1;
    for (i = 1; i <= N_max; i++)
    {        
        if (N_uis >= 4 && N_is < N_th)
        {
            // Variable for storing the number of stars matched in a particular iteration    
            int N_match = 0;

            // This will store the 4 stars used from the sm_3D_vecs table
            double four_stars[4][4];

            // Here the variable countt is used just to count whether 4 stars have been extracted
            int countt = 0, j = 0;

            // Pick body vectors of 4 stars from sm_3D_vecs
            for (countt = 0, j = 0; j < N_i && countt < 4; j++)
            {
                if ((int)sm_3D_vecs[j][0] != -1){
                    int k = 0;
                    for (k = 0; k < 4; k++)
                    {
                        four_stars[countt][k] = sm_3D_vecs[j][k];
                    }
                    countt++;
                }
            }
            int var=0;
            for(var = 0; var < 4; var++){
                printf("%lf  ",four_stars[var][0]);
            }
            // Perform 4 star N star algorithm
            sm_4_star(four_stars, sm_3D_vecs, sm_IS, body_vecs_IS, sm_K_vec_arr, &N_match, N_i, q, m, N_is);

            // Decrement matched stars from those detected from FE block for next iteration
            N_uis -= N_match;

            // Increment identified stars with matched stars
            N_is += N_match;

            // If no star is detected or we have not circulated the sm_3D_vecs exhaustively
            if (N_match == 0 && N_circ <= 2 * N_i){

                sm_4_star_circulate(sm_3D_vecs, &N_circ, N_i);
            
                if (N_circ >= 2 * N_i){
                    break;
                }
            }
        }

        // End the loop if the conditions "N_uis >= 4 && N_is < N_th" are not satisfied
        else{
            break;
        }
    }

    // sm_validate(sm_3D_vecs, sm_IS, body_vecs_IS, sm_GC, &N_is, N_i, tol, p_1, p_2);

    // Index variable for organising SM output
    int data_index = 0;

    double x_inertial[N_is], y_inertial[N_is], z_inertial[N_is];

    // Arrays to store body frame components (sm_3D_vec)
    double x_body[N_is], y_body[N_is], z_body[N_is];

    // If "some" stars match
    if (N_is != 0){
        *matched_stars = N_is;
    
        i = 0;
        for (i = 0; i < N_GC; i++){
            if (sm_IS[i][0] != -1){

                // Storing input ids and SSP ids in respective arrays
                input_ids[data_index] = sm_IS[i][0];
                star_ids[data_index] = sm_IS[i][1];

                // Storing body frame components of matched stars in respective arrays
                x_body[data_index] = body_vecs_IS[i][3];
                y_body[data_index] = body_vecs_IS[i][1];
                z_body[data_index] = body_vecs_IS[i][2];

                // Storinng intertial frame components of matched stars in respective arrays
                x_inertial[data_index] = sm_GC[sm_IS[i][1] - 1][1];
                y_inertial[data_index] = sm_GC[sm_IS[i][1] - 1][2];
                z_inertial[data_index] = sm_GC[sm_IS[i][1] - 1][3];

                // Increment data_index for next matched star
                data_index++;
            }
        }

        // Arranging SM output vectors compatible with Estimation block
        int j = 0;
        for (j = 0; j < N_is; j++){
            data[0][j] = x_inertial[j];
            data[0][j + N_is] = x_body[j];
            data[1][j] = y_inertial[j];
            data[1][j + N_is] = y_body[j];
            data[2][j] = z_inertial[j];
            data[2][j + N_is] = z_body[j];

        }
    }

    // If no star is matched
    else{
        printf("\nCOULD NOT MATCH STARS\n");
        
    }
}

//TM functions

void commonStars(int *common_stars, int prev_matched_stars, int curr_matched_stars, int prev_star_ids[prev_matched_stars], int prev_fe_ids[prev_matched_stars], int curr_star_ids[curr_matched_stars],
                 int curr_fe_ids[curr_matched_stars], int prev_tot_stars, int curr_tot_stars, double prev_centroids_st[MAX_STARS][3], double curr_centroids_st[MAX_STARS][3], double common_stars_data[MAX_STARS][5])
{
    int idx = 0;

    for(int i = 0; i < prev_matched_stars; i++)
    {
        for(int j = 0; j < curr_matched_stars; j++)
        {
            if(prev_star_ids[i] == curr_star_ids[j])
            {
                common_stars_data[idx][0] = curr_star_ids[j];
                for(int k = 0; k < prev_tot_stars; k++)
                {
                    if(prev_fe_ids[i] == prev_centroids_st[k][0])
                    {
                        common_stars_data[idx][1] = prev_centroids_st[k][1];
                        common_stars_data[idx][2] = prev_centroids_st[k][2];
                    }
                }

                for(int k = 0; k < curr_tot_stars; k++)
                {
                    if(curr_fe_ids[j] == curr_centroids_st[k][0])
                    {
                        common_stars_data[idx][3] = curr_centroids_st[k][1];
                        common_stars_data[idx][4] = curr_centroids_st[k][2];
                    }
                }
            *common_stars += 1;
            idx += 1;
            }
        }
    }
    return;
}

void predictCentroid(double common_centroid_data[][5], double predict_centroids_st[][3], int common_stars)
{
    double phi = 0, theta = 0, psi = 0;
    double u, v, f = FOCAL_LENGTH;
    double dun, dvn;
    double m[2*common_stars][3];

    for (int k = 0; k < common_stars; k++)
    {
        u = common_centroid_data[k][3];     //curr_centroid_x
        v = common_centroid_data[k][4];     //curr_centroid_y
        m[2*k][0] = u*v/f;
        m[2*k][1] = -f - u*u/f;
        m[2*k][2] = v;
        m[2*k+1][0] = f + v*v/f;
        m[2*k+1][1] = -u*v/f;
        m[2*k+1][2] = -u;
     }

    /*for (int k=0; k<2*N; k++)
    cout<<m[k][0]/f<<' '<<m[k][1]/f<<' '<<m[k][2]/f<<endl;*/

    double mTm[3][3] = {0,0,0,0,0,0,0,0,0};

    for (int k=0; k<2*common_stars; k++)
    {
        mTm[0][0] += m[k][0]*m[k][0];
        mTm[1][0] += m[k][1]*m[k][0];
        mTm[2][0] += m[k][2]*m[k][0];
        mTm[0][1] += m[k][0]*m[k][1];
        mTm[1][1] += m[k][1]*m[k][1];
        mTm[2][1] += m[k][2]*m[k][1];
        mTm[0][2] += m[k][0]*m[k][2];
        mTm[1][2] += m[k][1]*m[k][2];
        mTm[2][2] += m[k][2]*m[k][2];
    }

    /*for (int k=0; k<3; k++)
        cout<<mTm[k][0]/f/f<<' '<<mTm[k][1]/f/f<<' '<<mTm[k][2]/f/f<<endl;
    cout<<endl<<endl;*/

    double determinant = 0;
    double minv[3][3];      //(mTm)^c / |mTm|

    for(int i = 0; i < 3; i++)
    {
        determinant = determinant + (mTm[0][i]*(mTm[1][(i+1)%3]*mTm[2][(i+2)%3] - mTm[1][(i+2)%3]*mTm[2][(i+1)%3]));
    }

    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            minv[i][j] = ((mTm[(i+1)%3][(j+1)%3] * mTm[(i+2)%3][(j+2)%3]) - (mTm[(i+1)%3][(j+2)%3]*mTm[(i+2)%3][(j+1)%3]))/ determinant;
        }
    }

    /*for (int k=0; k<3; k++)
        cout<<minv[k][0]*f*f<<' '<<minv[k][1]*f*f<<' '<<minv[k][2]*f*f<<endl;
    cout<<endl<<endl;*/

    double m3[3][2*common_stars];
    for (int k = 0; k < 2*common_stars; k++)
    {        
        m3[0][k] = minv[0][0]*m[k][0] + minv[0][1]*m[k][1] + minv[0][2]*m[k][2];
        m3[1][k] = minv[1][0]*m[k][0] + minv[1][1]*m[k][1] + minv[1][2]*m[k][2];
        m3[2][k] = minv[2][0]*m[k][0] + minv[2][1]*m[k][1] + minv[2][2]*m[k][2];
    }

    /*for (int k=0; k<2*N; k++)
        cout<<m3[0][k]*f<<' '<<m3[1][k]*f<<' '<<m3[2][k]*f<<endl;
    cout<<endl<<endl;*/

    double col[2*common_stars];
    for (int k = 0; k < common_stars; k++)
    {
        col[2*k] = (common_centroid_data[k][3] - common_centroid_data[k][1]);
        col[2*k+1] = (common_centroid_data[k][4] - common_centroid_data[k][2]);
    }

    /*for (int k=0; k<2*N; k++)
        cout<<col[k]<<endl;
    cout<<endl<<endl;*/

    for (int k = 0; k < 2*common_stars; k++)
    {
        phi += m3[0][k]*col[k];
        theta += m3[1][k]*col[k];
        psi += m3[2][k]*col[k];
    }

    for (int k = 0; k < common_stars; k++)
    {
        u = common_centroid_data[k][1];
        v = common_centroid_data[k][2];

        dun = (u*v/f)*phi + (-f - u*u/f)*theta + v*psi;
        dvn = (f + v*v/f)*phi + (-u*v/f)*theta + (-u)*psi;

        predict_centroids_st[k][0] = common_centroid_data[k][0];        //star_id
        predict_centroids_st[k][1] = common_centroid_data[k][3] + dun;
        predict_centroids_st[k][2] = common_centroid_data[k][4] + dvn;
        
    }
}

// NOT USED YET
// sorting array of 'n' rows with respect to 'ind' column out of 2 coulums 
void bubbleSort2(double arr[][2], int n, int ind)
{
    int i, j;
    long double temp[2];
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n - i - 1; j++)
        {
            if (arr[j][ind] > arr[j + 1][ind])
            {
                // swap the elements
                temp[0] = arr[j][0];
                temp[1] = arr[j][1];

                arr[j][0] = arr[j + 1][0];
                arr[j][1] = arr[j + 1][1];

                arr[j + 1][0] = temp[0];
                arr[j + 1][1] = temp[1];
            }
        }
    }
}

// sorting array of 'n' rows with respect to 'ind' column out of 3 coulums 
void bubbleSort3(double arr[][3], int n, int ind)
{
    int i, j;
    long double temp[3];
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n - i - 1; j++)
        {
            if (arr[j][ind] > arr[j + 1][ind])
            {
                // swap the elements
                temp[0] = arr[j][0];
                temp[1] = arr[j][1];
                temp[2] = arr[j][2];

                arr[j][0] = arr[j + 1][0];
                arr[j][1] = arr[j + 1][1];
                arr[j][2] = arr[j + 1][2];

                arr[j + 1][0] = temp[0];
                arr[j + 1][1] = temp[1];
                arr[j + 1][2] = temp[2];
            }
        }
    }
}


double gc_id_angdist(int star_id_1, int star_id_2)
{
    double dot, norm1, norm2, id_dist;

    dot = sm_GC[star_id_1 - 1][1]*sm_GC[star_id_2 - 1][1] + sm_GC[star_id_1 - 1][2]*sm_GC[star_id_2 - 1][2] + sm_GC[star_id_1 - 1][3]*sm_GC[star_id_2 -1][3];
    norm1 = sqrt(sm_GC[star_id_1 - 1][1]*sm_GC[star_id_1 - 1][1] + sm_GC[star_id_1 - 1][2]*sm_GC[star_id_1 - 1][2] + sm_GC[star_id_1 - 1][3]*sm_GC[star_id_1 - 1][3]);
    norm2 = sqrt(sm_GC[star_id_2 - 1][1]*sm_GC[star_id_2 - 1][1] + sm_GC[star_id_2 - 1][2]*sm_GC[star_id_2 - 1][2] + sm_GC[star_id_2 - 1][3]*sm_GC[star_id_2 - 1][3]);

    id_dist = dot/(norm1*norm2);
    return id_dist;
}


double centroid_angdist(double x1, double y1, double x2, double y2)
{
    double dot, norm1, norm2, centroid_dist;
    double f = FOCAL_LENGTH;

    dot = x1*x2 + y1*y2 + f*f;
    norm1 = sqrt(x1*x1 + y1*y1 + f*f);
    norm2 = sqrt(x2*x2 + y2*y2 + f*f);

    centroid_dist = dot/(norm1*norm2);
    return centroid_dist;
}

int radiusBasedMatching(int common_stars, int next_tot_stars, double predicted_centroids_st[][3], double next_centroids_st[][3], double RBM_match[][4])
{   
    // bubbleSort3(predicted_centroids_st, common_stars, 0);
    // bubbleSort3(next_centroids_st, next_tot_stars, 0); 

    int sole_matched = 0;
    int matched = 0;
    int not_matched = 0;
    double r = RBM_RADIUS * pixel_size;
    // printf("r = %f\n", r);

    for (int i = 0; i < common_stars; i++)
    {
        int stars_in_range = 0;
        for (int j = 0; j < next_tot_stars; j++)
        {            
            double dx = absoluteValue(next_centroids_st[j][1] - predicted_centroids_st[i][1]);    
            double dy = absoluteValue(next_centroids_st[j][2] - predicted_centroids_st[i][2]);

            // printf("x1 = %f, x2 = %f, dx = %f\n", next_centroids_st[j][1], predicted_centroids_st[i][1], dx);
            // double distance = sqrt(dx*dx + dy*dy);
            // printf("distance %f\n", distance);

            if (dx < r)
            {
                if (dy < r)
                {
                    stars_in_range++;
                    if (stars_in_range > 1)
                    {
                        //store not matched stars BUT NOT ALL ENTRIES ARE CONSIDERED IN THIS
                        // newEntries[not_matched][0] = predicted_centroids_st[i][0]; 
                        // newEntries[not_matched][1] = predicted_centroids_st[i][1]; 
                        // newEntries[not_matched][2] = predicted_centroids_st[i][2];
                        not_matched++; 
                        break;
                    }
                    sole_matched = j;
                }                                                    
            }
        }
        if (stars_in_range == 1)
        {
            // printf("check: stars_in_range = 1\n");
            RBM_match[matched][0] = next_centroids_st[sole_matched][0]; //fe_id
            RBM_match[matched][1] = predicted_centroids_st[i][0];       //star_id
            RBM_match[matched][2] = next_centroids_st[sole_matched][1]; //next_cent_x
            RBM_match[matched][3] = next_centroids_st[sole_matched][2]; //next_cent_y
            matched++;        
        }
    }
    return (matched);
}

void starNeighbourhoodMatch(double RBM_matched[][4], int RBM_matched_stars, double next_centroids_st[][3],int next_tot_stars,int sm_SNT[N_GC][17],double sm_GC[][4],double newEntries[][4], int* new_matched_stars)
{
    double fe_unmatched[next_tot_stars][3];
    int num_unmatched = 0;

    for (int j = 0; j < next_tot_stars; j++)
    {
        bool already_matched = false;

        for (int i = 0; i < next_tot_stars; i++)
        {
            if((RBM_matched[i][2] == next_centroids_st[j][1]) && (RBM_matched[i][3] == next_centroids_st[j][2]))
            {
                already_matched = true;
                break;
            }
        }
        if(!(already_matched))
        {
            fe_unmatched[num_unmatched][0] = next_centroids_st[j][0];
            fe_unmatched[num_unmatched][1] = next_centroids_st[j][1];
            fe_unmatched[num_unmatched][2] = next_centroids_st[j][2];
            // printf("%f %f %f\n", fe_unmatched[num_unmatched][0],fe_unmatched[num_unmatched][1],fe_unmatched[num_unmatched][2]);
            num_unmatched++;
        }
    }

    int i_unmatch = 0;
    int matched_stars = 0;

    for (i_unmatch = 0; i_unmatch < num_unmatched; i_unmatch++)
    {
        bool done = false;
        for (int i = 0; i < RBM_matched_stars; i++)
        {
            int curr_ref_star = RBM_matched[i][1];  //star_id
            // printf("curr_star_id%d\n", curr_ref_star);

            if (curr_ref_star > 5060) continue;     //becz we don't have SNT for N_GC=8876 yet
            if (matched_stars == Nth - next_tot_stars) break;

            int j = 0;
            while(sm_SNT[curr_ref_star - 1][j] != 0)
            {
                double x1 = RBM_matched[i][2];
                double y1 = RBM_matched[i][3];
                double x2 = fe_unmatched[i_unmatch][1];
                double y2 = fe_unmatched[i_unmatch][2];

                double GC_ID_angdist = gc_id_angdist(curr_ref_star, sm_SNT[curr_ref_star - 1][j]);
                double Centroid_angdist = centroid_angdist(x1, y1, x2, y2);
                double error = absoluteValue(GC_ID_angdist - Centroid_angdist);

                // if (10e6 * error < 100)
                // printf("GC_angdist =%f, Centroid_angdist =%f, Error =%f\n", GC_ID_angdist, Centroid_angdist, 10e6 * error);

                if(10e6*error < 2)
                {
                    newEntries[matched_stars][0] = fe_unmatched[i_unmatch][0];
                    newEntries[matched_stars][1] = sm_SNT[curr_ref_star - 1][j];
                    newEntries[matched_stars][2] = fe_unmatched[i_unmatch][1];
                    newEntries[matched_stars][3] = fe_unmatched[i_unmatch][2];
                    // printf("matched_stars =%d SNT_star = %d ", matched_stars, sm_SNT[curr_ref_star - 1][j]);
                    // printf("fe_id =%f, star_id =%f \n\n", newEntries[matched_stars][0], newEntries[matched_stars][1]);
                    matched_stars++;
                    done = true;
                }
                j++;
                if (done) break;

            }
            if (done) break;
        }
    }
    *new_matched_stars = matched_stars;
}