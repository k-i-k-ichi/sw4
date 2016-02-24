#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define EPSI 0.0000000001

typedef struct vertice{
    double x;
    double y;
    int isNull;
}Vertice;


void set_null(Vertice **a)
{
    *a = NULL;
}

// Edges of polygon are a linked list of pair of vertices
typedef struct edge{
    struct edge* prevEdge;
    struct edge* nextEdge;
    Vertice* a;
    Vertice* b;
}
    
    
   


double polar_angle(double p0_x, double p0_y, double p1_x, double p1_y){
    double deltaY = p2_y - p1_y;
    double deltaX = p2_x - p1_x;
    double polar_angle = atan2(deltaY/deltaX);
    return polar_angle;
}



// for qsort()
int md_comparator(const void *v1, const void *v2)
{
    const Vertice* p1 = (Vertice* )v1;
    const Vertice* p2 = (Vertice* )v2;
    if ((p1->angle_to_guard - p2->angle_to_guard) < EPSI) // epsilon equality
        return -1;
    else if ((p1->angle_to_guard - p2->angle_to_guard) > EPSI)
        return +1;
    else
        return 0;
}




int get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y, double* i_x, double* i_y)
{
    double s1_x, s1_y, s2_x, s2_y;
    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

    double s, t;
    s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
    t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
    {
        // Collision detected
        if (i_x != NULL)
            *i_x = p0_x + (t * s1_x);
        if (i_y != NULL)
            *i_y = p0_y + (t * s1_y);
        return 1;
    }
    return 0; // No collision
}




int is_inside(Vertice* guard, Vertice* polygon){
    // cast a ray vertically downward from guard
    // as a line segment from guard to a very large negative vertex 
    Vertice* temp = malloc(sizeof(struct vertice));
    temp->x = guard->x;
    temp->y = -999999;
    
    int intersectCount =0;
    
    int vertexIndex = 0;
    while(polygon[vertexIndex+1].isNull ==0){ 
        // If guard stand directly on top of vertex
        if(guard->x -
    
        // Iterate through all of edges of polygon
        // Check for collision with the downward ray
        int isIntersect = get_line_intersection(guard->x, guard->y, temp->x, temp->y, polygon[vertexIndex].x, polygon[vertexIndex].y, polygon[vertexIndex+1].x, polygon[vertexIndex+1].y, NULL, NULL);
        
        if (isIntersect == 1) intersectCount++;
        vertexIndex++;
    }
    
    // the last edge is from last vertex to first vertex
    int isIntersect = get_line_intersection(guard->x, guard->y, temp->x, temp->y, polygon[vertexIndex].x, polygon[vertexIndex].y, polygon[0].x, polygon[0].y, NULL, NULL);
    if(isIntersect==1)intersectCount++;
}

 

Vertice* parse(char* buffer){
    int memberCount=0;
    int i=0;
    while(buffer[i]!='\0'){
        if(buffer[i] == '(') memberCount++;
        i++;
    }
    Vertice* verticeArray = malloc( (memberCount+1) *sizeof(struct vertice));
    
    int pointer=0;
    for(i=0; i< memberCount;i++){
        
        char doubleBuffer[100];
        while( !isdigit(buffer[pointer]) && buffer[pointer] != '-') pointer++; // seek to first digit or minus sign
        
        
        int initialPointer = pointer;
        while(buffer[pointer]!= ','){ // x coordinate
            doubleBuffer[pointer-initialPointer] = buffer[pointer];
            pointer++;
        }
        
        while( !isdigit(buffer[pointer]) && buffer[pointer] != '-') pointer++; // seek to next digit or minus sign
        
        
        char doubleBuffer2[100];
        initialPointer = pointer;
        while(buffer[pointer]!=')'){ // y coordinate
            doubleBuffer2[pointer-initialPointer] = buffer[pointer];
            pointer++;
        }
        
        
        
        Vertice* new = malloc(sizeof(struct vertice));
        new->x = atof(doubleBuffer);
        new->y = atof(doubleBuffer2);
        new->isNull =0;
        verticeArray[i] = *new;
    }
    
    free(buffer);
    verticeArray[i+1].x=0, verticeArray[i+1].y=0, verticeArray[i+1].isNull=1;
    
    return verticeArray;
}


int main(){
    FILE* fp;
    fp = fopen("check.pol", "r");
    for(int main=0; main<1; main++){
        char* buffer = malloc(25000*sizeof(char));
        fgets(buffer, 25000, fp);
        
        puts(buffer);
        
        int i =0;
        while(buffer[i] != '(') i++; // seek first parenth        
        
        char* polygon = malloc(20000*sizeof(char));
        
        int starter = i;
        while(buffer[i] != ';'){
             polygon[i-starter] = buffer[i];
             i++;
        }
        polygon[i-starter]='\0';
        
        
        puts(polygon);
        
        while( buffer[i] != '(')i++; // seek to next parenth
        
        starter = i;
        char* guard = malloc(8000*sizeof(char));
        while(buffer[i] != '\n'){
                guard[i-starter] = buffer[i];
                i++;
        }
        guard[i-starter]='\0';
        puts(guard);
        free(buffer);
        
        Vertice* polygonArray = parse(polygon);
        Vertice* guardArray = parse(guard);
         
    } 
}
