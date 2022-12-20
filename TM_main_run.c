#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

//images
#include "image_1.h"
#include "image_2.h" 
#include "image_3.h" 
#include "image_4.h" 
#include "image_5.h" 
#include "image_6.h" 
#include "image_7.h" 
#include "image_8.h" 
#include "image_9.h" 
#include "image_10.h"

#include "functions.h"
#include "constants.h"

short (*frames[100]) [810] = 
{arr_out_img_1, arr_out_img_2, arr_out_img_3, arr_out_img_4, arr_out_img_5, arr_out_img_6, arr_out_img_7, arr_out_img_8, arr_out_img_9, arr_out_img_10,};

int sm_K_vec_arr[224792][3];

int main()
{
    int frame = 1;
    bool TM_on = false;

    //select previous frame
        int prev_frame = frame - 1;
        short (*prev_img) [810] = frames[prev_frame];
        int prev_tot_stars = 0; 
        double prev_centroids_st[MAX_STARS][3];    
        
    //FE on previous frame
        FE(prev_img, prev_centroids_st, &prev_tot_stars);

        int prev_matched_stars = 0;
        int prev_input_ids[MAX_STARS] = {0};
        int prev_star_ids[MAX_STARS] = {0};
        double prev_data[3][MAX_STARS];

    //LISM on previous frame
        LISM(prev_centroids_st, prev_tot_stars, prev_data, prev_input_ids, prev_star_ids, &prev_matched_stars);

    //---------------------------------------------------------------------------------------------------------------------------
    //select current frame
        int curr_frame = frame;
        short (*curr_img) [810] = frames[curr_frame];
        int curr_tot_stars = 0;
        double curr_centroids_st[MAX_STARS][3];        

    //FE on the current frame
        FE(curr_img, curr_centroids_st, &curr_tot_stars);

        int curr_matched_stars = 0;
        int curr_input_ids[MAX_STARS] = {0};
        int curr_star_ids[MAX_STARS] = {0};
        double curr_data[3][MAX_STARS];

    //LISM on current frame
        LISM(curr_centroids_st, curr_tot_stars, curr_data, curr_input_ids, curr_star_ids, &curr_matched_stars);

    //----------------------------------------------------------------------------------------------------------------------------

    //select next frame
        int next_frame = frame - 1;
        int next_tot_stars = 0; 
        double next_centroids_st[MAX_STARS][3];

        int next_matched_stars = 0;
        int next_input_ids[MAX_STARS] = {0};
        int next_star_ids[MAX_STARS] = {0};
        double next_data[3][MAX_STARS]; 

    while(frame < 10){

        if (frame > 1)
        {
            //pointing prev frame to curr frame
            prev_tot_stars = curr_tot_stars;
            prev_matched_stars = curr_matched_stars;

            for(int i = 0; i < MAX_STARS; i++)
            {
                //here both fe_id, star_id are of same ith star becz these are assigned in LISM
                prev_input_ids[i] = curr_input_ids[i];
                prev_star_ids[i] = curr_star_ids[i];

                //this ith centroids need not belong to ith star
                prev_centroids_st[i][0] = curr_centroids_st[i][0];
                prev_centroids_st[i][1] = curr_centroids_st[i][1];
                prev_centroids_st[i][2] = curr_centroids_st[i][2];

                //curr_data -> prev_data 
                for (int i = 0; i < MAX_STARS; i++)
                {
                    prev_data[0][i] = curr_data[0][i];
                    prev_data[1][i] = curr_data[1][i];
                    prev_data[2][i] = curr_data[2][i];
                }               
            }

            //fe_centroids next_frame -> curr_frame
            curr_tot_stars = next_tot_stars;

            for(int i = 0; i < MAX_STARS; i++)
            {
                curr_centroids_st[i][0] = next_centroids_st[i][0];
                curr_centroids_st[i][1] = next_centroids_st[i][1];
                curr_centroids_st[i][2] = next_centroids_st[i][2];
            }

            if(! TM_on)
            {   
                printf("TM Failed\n");
                curr_matched_stars = 0;
                
                //LISM on current frame
                LISM(curr_centroids_st, curr_tot_stars, curr_data, curr_input_ids, curr_star_ids, &curr_matched_stars);
            }
            else
            {
                curr_matched_stars = next_matched_stars;

                for(int i = 0; i < MAX_STARS; i++)
                {
                    curr_input_ids[i] = next_input_ids[i];
                    curr_star_ids[i] = next_star_ids[i];
                    //next_data -> curr_data
                    for (int i = 0; i < MAX_STARS; i++)
                    {
                        curr_data[0][i] = next_data[0][i];
                        curr_data[1][i] = next_data[1][i];
                        curr_data[2][i] = next_data[2][i];
                    }
                }                
            }     
        }

        //variable to store common stars in both the frames
        int common_stars = 0;

        //array to store commmon stars
        double common_centroid_data[MAX_STARS][5];

        commonStars(&common_stars, prev_matched_stars, curr_matched_stars, prev_star_ids, prev_input_ids, curr_star_ids, curr_input_ids, prev_tot_stars, curr_tot_stars, prev_centroids_st, curr_centroids_st, common_centroid_data);

        printf("Common Stars: %d\n", common_stars);

        //--------------------------------------------------------------------------------------------------------------------

        double predicted_centroids_st[common_stars][3];

        if(common_stars >= 2)
        {
            //predict next frame centroids based on the data of previous frame and next frame
            predictCentroid(common_centroid_data, predicted_centroids_st, common_stars);

            next_frame = frame + 1;
            frame ++;
            next_tot_stars = 0;
            short (*next_img) [810] = frames[next_frame];
                       
            //FE on the next frame
            FE(next_img, next_centroids_st, &next_tot_stars);

            double RBM_centroid_st[MAX_STARS][4];
            int RBM_matched_stars = radiusBasedMatching(common_stars, next_tot_stars, predicted_centroids_st, next_centroids_st, RBM_centroid_st);
            printf("Matched stars = %d ", RBM_matched_stars);

            //next fe, star IDs of matched stars
            for (int i = 0; i < RBM_matched_stars; i++)
            {
                next_input_ids[i] = RBM_centroid_st[i][0];
                next_star_ids[i] = RBM_centroid_st[i][1];
            }

            if(RBM_matched_stars > N_TH_TRACKING)
            {  
                //next data from sm_GC
                for (int i = 0; i < MAX_STARS; i++)
                {
                    next_data[0][i] = sm_GC[next_star_ids[i]][1];
                    next_data[1][i] = sm_GC[next_star_ids[i]][2];
                    next_data[2][i] = sm_GC[next_star_ids[i]][3];
                }

                next_matched_stars = RBM_matched_stars;
                TM_on = true;
                printf("GO TO ESTIMAION\n");
            }

            else{

                //need to update next_matched_stars, next_data
                // from newEntries & RBM_match
                int new_matched_stars = 0;
                double newEntries[N_TH_TRACKING][3];
                starNeighbourhoodMatch(RBM_centroid_st, next_centroids_st, next_tot_stars, sm_SNT, sm_GC, newEntries, &new_matched_stars);
                
                next_matched_stars = RBM_matched_stars + new_matched_stars;
                printf("Identify new stars entering the FOV\n");
            }

        }
        else{
            frame ++;
            TM_on = false;
            continue;
        }

    }    
}