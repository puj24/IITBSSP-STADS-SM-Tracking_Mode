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
#include "image_11.h"
#include "image_12.h"
#include "image_13.h"
#include "image_14.h"
#include "image_15.h"
#include "image_16.h"
#include "image_17.h"
#include "image_18.h"
#include "image_19.h"
#include "image_20.h"
#include "image_21.h"
#include "image_22.h"
#include "image_23.h"
#include "image_24.h"
#include "image_25.h"
#include "image_26.h"
#include "image_27.h"
#include "image_28.h"
#include "image_29.h"
#include "image_30.h"
#include "image_31.h"
#include "image_32.h"
#include "image_33.h"
#include "image_34.h"
#include "image_35.h"
#include "image_36.h"
#include "image_37.h"
#include "image_38.h"
#include "image_39.h"
#include "image_40.h"
#include "image_41.h"
#include "image_42.h"
#include "image_43.h"
#include "image_44.h"
#include "image_45.h"
#include "image_46.h"
#include "image_47.h"
#include "image_48.h"
#include "image_49.h"
#include "image_50.h"
#include "image_51.h"
#include "image_52.h"
#include "image_53.h"
#include "image_54.h"
#include "image_55.h"
#include "image_56.h"
#include "image_57.h"
#include "image_58.h"
#include "image_59.h"
#include "image_60.h"
#include "image_61.h"
#include "image_62.h"
#include "image_63.h"
#include "image_64.h"
#include "image_65.h"
#include "image_66.h"
#include "image_67.h"
#include "image_68.h"
#include "image_69.h"
#include "image_70.h"
#include "image_71.h"
#include "image_72.h"
#include "image_73.h"
#include "image_74.h"
#include "image_75.h"
#include "image_76.h"
#include "image_77.h"
#include "image_78.h"
#include "image_79.h"
#include "image_80.h"
#include "image_81.h"
#include "image_82.h"
#include "image_83.h"
#include "image_84.h"
#include "image_85.h"
#include "image_86.h"
#include "image_87.h"
#include "image_88.h"
#include "image_89.h"
#include "image_90.h"
#include "image_91.h"
#include "image_92.h"
#include "image_93.h"
#include "image_94.h"
#include "image_95.h"
#include "image_96.h"
#include "image_97.h"
#include "image_98.h"
#include "image_99.h"
#include "image_100.h"

#include "functions.h"
#include "constants.h"

short (*frames[100]) [810] = 
{arr_out_img_1, arr_out_img_2, arr_out_img_3, arr_out_img_4, arr_out_img_5, arr_out_img_6, arr_out_img_7, arr_out_img_8, arr_out_img_9, arr_out_img_10,
 arr_out_img_11, arr_out_img_12, arr_out_img_13, arr_out_img_14, arr_out_img_15, arr_out_img_16, arr_out_img_17, arr_out_img_18, arr_out_img_19, arr_out_img_20,
 arr_out_img_21, arr_out_img_22, arr_out_img_23, arr_out_img_24, arr_out_img_25, arr_out_img_26, arr_out_img_27, arr_out_img_28, arr_out_img_29, arr_out_img_30,
 arr_out_img_31, arr_out_img_32, arr_out_img_33, arr_out_img_34, arr_out_img_35, arr_out_img_36, arr_out_img_37, arr_out_img_38, arr_out_img_39, arr_out_img_40,
 arr_out_img_41, arr_out_img_42, arr_out_img_43, arr_out_img_44, arr_out_img_45, arr_out_img_46, arr_out_img_47, arr_out_img_48, arr_out_img_49, arr_out_img_50, 
 arr_out_img_51, arr_out_img_52, arr_out_img_53, arr_out_img_54, arr_out_img_55, arr_out_img_56, arr_out_img_57, arr_out_img_58, arr_out_img_59, arr_out_img_60,
 arr_out_img_61, arr_out_img_62, arr_out_img_63, arr_out_img_64, arr_out_img_65, arr_out_img_66, arr_out_img_67, arr_out_img_68, arr_out_img_69, arr_out_img_70,
 arr_out_img_71, arr_out_img_72, arr_out_img_73, arr_out_img_74, arr_out_img_75, arr_out_img_76, arr_out_img_77, arr_out_img_78, arr_out_img_79, arr_out_img_80, 
 arr_out_img_81, arr_out_img_82, arr_out_img_83, arr_out_img_84, arr_out_img_85, arr_out_img_86, arr_out_img_87, arr_out_img_88, arr_out_img_89, arr_out_img_90, 
 arr_out_img_91, arr_out_img_92, arr_out_img_93, arr_out_img_94, arr_out_img_95, arr_out_img_96, arr_out_img_97, arr_out_img_98, arr_out_img_99, arr_out_img_100};

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
        printf("Frame1: calling LISM\n\n");
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
        printf("Frame2: calling LISM\n\n");
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

    while(frame < 100){

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
                printf("Frame%d: TM Failed. Calling LISM\n\n", frame + 1);
                curr_matched_stars = 0;
                
                //LISM on current frame
                LISM(curr_centroids_st, curr_tot_stars, curr_data, curr_input_ids, curr_star_ids, &curr_matched_stars);
            }
            else
            {   
                printf("Frame%d: TM succees\n\n", frame + 1);
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

        if (frame == 99) break;
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

            // for(int i = 0; i <  common_stars; i++)
            // {   
            //         printf("%.16f %.16f %.16f\n",  predicted_centroids_st[i][0], predicted_centroids_st[i][1], predicted_centroids_st[i][2]);                
            // }

            next_frame = frame + 1;
            frame ++;
            next_tot_stars = 0;
            short (*next_img) [810] = frames[next_frame];
                       
            //FE on the next frame
            FE(next_img, next_centroids_st, &next_tot_stars);

            // for(int i = 0; i < next_tot_stars; i++){
            //     printf("%d %.16f %.16f\n", (int)next_centroids_st[i][0], next_centroids_st[i][1], next_centroids_st[i][2]);
            // }

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
                double newEntries[N_TH_TRACKING][4];
                starNeighbourhoodMatch(RBM_centroid_st, next_centroids_st, next_tot_stars, sm_SNT, sm_GC, newEntries, &new_matched_stars);
                
                printf("newMatched_stars: %d\n", new_matched_stars);
                next_matched_stars = RBM_matched_stars + new_matched_stars;

                for (int i = RBM_matched_stars; i < next_matched_stars; i++)
                {
                    next_input_ids[i] = newEntries[i - RBM_matched_stars][0];
                    next_star_ids[i] = newEntries[i - RBM_matched_stars][1];
                }

                for (int i = 0; i < MAX_STARS; i++)
                {
                    next_data[0][i] = sm_GC[next_star_ids[i]][1];
                    next_data[1][i] = sm_GC[next_star_ids[i]][2];
                    next_data[2][i] = sm_GC[next_star_ids[i]][3];
                }
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