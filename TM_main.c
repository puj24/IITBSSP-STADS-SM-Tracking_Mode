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
    while(frame < 100){

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

        //-----------------------------------------------------------------------------------------------------------------

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

        //--------------------------------------------------------------------------------------------------------------------

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

            int next_frame = frame + 1;
            short (*next_img) [810] = frames[next_frame];
            int next_tot_stars = 0;
            double next_centroids_st[MAX_STARS][3];

            //FE on the next frame
            FE(next_img, next_centroids_st, &next_tot_stars);

            double RBM_centroid_st[next_tot_stars][4];
            int matched_stars = radiusBasedMatching(common_stars, next_tot_stars, predicted_centroids_st, next_centroids_st, RBM_centroid_st);
            printf("Matched stars = %d ", matched_stars);

            if(matched_stars > N_TH_TRACKING)
            {  
                
                printf("GO TO ESTIMAION");
            }

            else{
                int new_matched_stars = 0;
                double newEntries[N_TH_TRACKING][4];
                starNeighbourhoodMatch(RBM_centroid_st, next_centroids_st, next_tot_stars, sm_SNT, sm_GC, newEntries, &new_matched_stars);
                printf("Identify new stars entering the FOV");
            }

        }
        else{
            continue;
        }

    }    
}