#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define a_1 1
#define a_2 2

int main(){
    int p;

    // for(int i = 1; i < 3; i++)
    // {
    //     p = dprintf("a_%d", i);
    // }
    // printf("%d", p);
    // printf(sprintf("a_%d", 1));

    for (int i = 0; i < 100; i++)
    {
        printf("tot_stars_UIS[%d] = N_i_%d;\n", i, i + 1);
        // tot_stars_UIS[0] = N_i_1;
    }
    
}