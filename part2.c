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
    double angle_to_guard;
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
}Edge;

double polar_angle(double p1_x, double p1_y, double p2_x, double p2_y){
    double deltaY = p2_y - p1_y;
    double deltaX = p2_x - p1_x;
    double polar_angle = atan2(deltaY, deltaX);
    return polar_angle;
}



// for qsort()
int md_comparator(const void *v1, const void *v2)
{
    const Vertice* p1 = (Vertice* )v1;
    const Vertice* p2 = (Vertice* )v2;
    if ((p1->angle_to_guard - p2->angle_to_guard) < -EPSI) // epsilon equality
        return -1;
    else if ((p1->angle_to_guard - p2->angle_to_guard) > EPSI)
        return +1;
    else
        return 0;
}




int get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y, double* i_x, double* i_y)
{

    // special case
    double case1, case2, case3, case4;
    case1 = fabs(p0_x-p2_x) + fabs(p0_y-p2_y);
    case2 = fabs(p0_x-p3_x) + fabs(p0_y-p3_y);
    case3 = fabs(p1_x-p2_x) + fabs(p1_y-p2_y);
    case4 = fabs(p1_x-p3_x) + fabs(p1_y-p3_y);
    
    if(case1<2*EPSI || case2<2*EPSI){
        if(i_x != NULL)*i_x = p0_x;
        return 1;
    }
    
    if(case3<2*EPSI || case4<2*EPSI){
        if(i_x != NULL)*i_x = p1_x;
        return 1;
    }
    
        
    
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
    //paraemter: guard : pointer to guard
    //      polygon: pointer to array of polygon
    
    // cast a ray vertically downward from guard
    // as a line segment from guard to a very large negative vertex 
    Vertice* temp = malloc(sizeof(struct vertice));
    temp->x = guard->x;
    temp->y = -999999;
    
    int intersectCount =0;
    
    int vertexIndex = 0;
    while(polygon[vertexIndex+1].isNull ==0){
        
        if (fabs(guard->x - polygon[vertexIndex].x) <EPSI && fabs(guard->y - polygon[vertexIndex].y) < EPSI ) return -1;
        // Iterate through all of edges of polygon
        // Check for collision with the downward ray
        int isIntersect = get_line_intersection(guard->x, guard->y, temp->x, temp->y, polygon[vertexIndex].x, polygon[vertexIndex].y, polygon[vertexIndex+1].x, polygon[vertexIndex+1].y, NULL, NULL);
        
        if (isIntersect == 1) intersectCount++;
        vertexIndex++;
    }
    
    // the last edge is from last vertex to first vertex
    

    if (fabs(guard->x - polygon[vertexIndex].x) <EPSI && fabs(guard->y - polygon[vertexIndex].y) < EPSI ) return -1;
    int isIntersect = get_line_intersection(guard->x, guard->y, temp->x, temp->y, polygon[vertexIndex].x, polygon[vertexIndex].y, polygon[0].x, polygon[0].y, NULL, NULL);
    if(isIntersect==1)intersectCount++;
    
    free(temp);
    
    if(intersectCount %2 ==0) return 0; // even intersect
    else return 1;
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
        
        char doubleBuffer[100]={0};
        while( !isdigit(buffer[pointer]) && buffer[pointer] != '-') pointer++; // seek to first digit or minus sign
        
        
        int initialPointer = pointer;
        while(buffer[pointer]!= ','){ // x coordinate
            doubleBuffer[pointer-initialPointer] = buffer[pointer];
            pointer++;
        }
        
        while( !isdigit(buffer[pointer]) && buffer[pointer] != '-') pointer++; // seek to next digit or minus sign
        
        
        char doubleBuffer2[100]={0};
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
    verticeArray[memberCount].x=0, verticeArray[memberCount].y=0, verticeArray[memberCount].isNull=1; // mark the end with a null vertice
    
    return verticeArray;
}


int main(){
    FILE* fp;
    fp = fopen("check.pol", "r");
    for(int main=0; main<20; main++){
        char* buffer = calloc(30000,sizeof(char));
        fgets(buffer, 30000, fp);
        
        
        int i =0;
        while(buffer[i] != '(') i++; // seek first parenth        
        
        char* polygon = calloc(20000,sizeof(char));
        
        int starter = i;
        while(buffer[i] != ';'){
             polygon[i-starter] = buffer[i];
             i++;
        }
        polygon[i-starter]='\0';
        
        
        
        
        while( buffer[i] != '(')i++; // seek to next parent(
        
        starter = i;
        char* guard = calloc(10000,sizeof(char));
        while(buffer[i] != '\n' && buffer[i] != '\0'){
                guard[i-starter] = buffer[i];
                i++;
        }
        guard[i-starter]='\0';
        
        free(buffer);
        
        Vertice* polygonArray = parse(polygon);
        Vertice* guardArray = parse(guard);
        int k=0;
        printf("main: %d\n", main);
        while (guardArray[k].isNull==0){
            printf("guard %d : inside: %d\n",k, is_inside(&guardArray[k], polygonArray));
            k++;
        }
        double a, b;
        free(polygonArray);
        free(guardArray);
    } 
}
