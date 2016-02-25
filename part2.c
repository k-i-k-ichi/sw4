#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define EPSI 0.0000000001
#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)
#define INSIDE 0
#define OUTSIDE 1

typedef struct vertice{
    double x;
    double y;
    int isNull;
    double angle_to_guard;
    int isInside;
    int isVisible;
}Vertice;


void set_null(Vertice **a)
{
    *a = NULL;
}

// Edges of polygon are a linked list of pair of vertices
typedef struct edge{
    struct edge* prevEdge;
    struct edge* nextEdge;
    Vertice* strt;
    Vertice* fins;
    int check;
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
    
    // sort by angle 
    if ((p1->angle_to_guard - p2->angle_to_guard) < -EPSI) // epsilon equality
        return -1;
    else if ((p1->angle_to_guard - p2->angle_to_guard) > EPSI)
        return +1;
    
    // angle are equal, sort by order
    if(p1->isNull > p2->isNull) // Greater isNull means vertex p1 is before p2 in poylgon array 
        return -1;
    else if (p1->isNull < p2->isNull)
        return 1;
    else // 2 vertice have the same order, must be the same vetex
        return 0;

}




int get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y, double* i_x, double* i_y)
{

    // Return: 0 for non collision
    //          1 for proper crossing
    //            -1 when they are "touching"

    // special case
    double case1, case2, case3, case4;
    case1 = fabs(p0_x-p2_x) + fabs(p0_y-p2_y);
    case2 = fabs(p0_x-p3_x) + fabs(p0_y-p3_y);
    case3 = fabs(p1_x-p2_x) + fabs(p1_y-p2_y);
    case4 = fabs(p1_x-p3_x) + fabs(p1_y-p3_y);
    
    if(case1<2*EPSI || case2<2*EPSI){
        if(i_x != NULL)*i_x = p0_x;
        if(i_y != NULL)*i_y = p0_y;
        return -1;
    }
    
    if(case3<2*EPSI || case4<2*EPSI){
        if(i_x != NULL)*i_x = p1_x;
        if(i_y != NULL)*i_y = p1_y;
        return -1;
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
        double temp_x = p0_x + (t * s1_x);
        double temp_y = p0_y + (t * s1_y);
        if (i_x != NULL)    
            *i_x = temp_x;
        if (i_y != NULL)
            *i_y = temp_y;
        if (fabs(temp_x-p0_x)< EPSI && fabs(temp_y-p0_y)<EPSI) return -1;
        if (fabs(temp_x-p1_x)< EPSI && fabs(temp_y-p1_y)<EPSI) return -1;
        if (fabs(temp_x-p2_x)< EPSI && fabs(temp_y-p2_y)<EPSI) return -1;
        if (fabs(temp_x-p3_x)< EPSI && fabs(temp_y-p3_y)<EPSI) return -1;
        return 1;
    }
    return 0; // No collision, actually never happens
}




int is_inside(Vertice p, Vertice* polygon){
    int counter = 0;
    int i;
    double xinters;
    Vertice p1,p2;
    int N = polygon[0].isNull;
    p1 = polygon[0];
    for (i=1;i<=N;i++) {
        if (fabs(p.x-p1.x)< EPSI && fabs(p.y-p1.y) <EPSI) return -1;
        
        p2 = polygon[i % N];
        if (p.y > MIN(p1.y,p2.y)-EPSI) {
            if (p.y <= MAX(p1.y,p2.y)+EPSI) {
                if (p.x <= MAX(p1.x,p2.x)+EPSI) {
                    if (fabs(p1.y - p2.y)<EPSI) {
                        xinters = (p.y-p1.y)*(p2.x-p1.x)/(p2.y-p1.y)+p1.x;
                        if (p1.x == p2.x || p.x <= xinters) counter++;
                        
                    }
                }
            }
        }
    
        p1 = p2;
    }

    if (counter % 2 == 0)
        return(OUTSIDE);
    else
        return(INSIDE);
}
        

 

Vertice* parse(char* buffer){
    int memberCount=0;
    int i=0;
    while(buffer[i]!='\0'){
        if(buffer[i] == '(') memberCount++;
        i++;
    }
    Vertice* verticeArray = calloc((memberCount+1),sizeof(struct vertice));
    
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
        new->isNull =memberCount-i;
        new->isVisible = 0;
        verticeArray[i] = *new;
    }
    
    free(buffer);
    verticeArray[memberCount].x=0, verticeArray[memberCount].y=0, verticeArray[memberCount].isNull=0; // mark the end with a null vertice
    verticeArray[memberCount].isVisible = 0;
    
    return verticeArray;
}

void copy_vertice_array(Vertice* dest, Vertice* source){
    for (int i =0; i<source[0].isNull; i++){
        Vertice* new = malloc(sizeof(struct vertice));
        new->x = source[i].x; new->y = source[i].y;
        new->isNull = source[i].isNull; new->angle_to_guard = source[i].angle_to_guard;
        new->isInside = source[i].isInside;
        new->isVisible = source[i].isVisible;
        dest[i] = *new;
    }
}


int main(){
    FILE* fp;
    fp = fopen("check.pol", "r");
    for(int main=0; main<1; main++){
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

        // List all guard inside polygon
        int k=0;
        int counter = 0;
        while (guardArray[k].isNull!=0){
            
            if(is_inside(guardArray[k], polygonArray) != 0){
                guardArray[k].isInside =1;
                counter++;
            }
            k++;
            
        }
        Vertice* insideGuardArray = calloc(counter,sizeof(struct vertice));
        i=0;
        for(k = 0; k<guardArray[0].isNull; k++){
            if(guardArray[k].isInside==1){
                insideGuardArray[i]=guardArray[k];
                insideGuardArray[i].isNull=counter-i;
                i++;
            }
        }

        
        // Create edge linked list
        Edge* first = malloc(sizeof(struct edge));
        first->prevEdge = NULL; first->nextEdge = NULL; first->strt = &polygonArray[0]; first->fins = &polygonArray[1]; first->check =0;
        
        // for each 2 vertices of polygon make an edge object
        Edge* edgePointer = first;
        for (int k =1; k<polygonArray[0].isNull-1; k++) {
            printf(" ! %d\n",k);
            Edge* new = malloc(sizeof(struct edge));
            new->strt = &polygonArray[k];
            new->fins = &polygonArray[k+1];
            new->prevEdge = edgePointer;
            new->check=0;
            edgePointer->nextEdge = new;
            edgePointer = new;
        }
        edgePointer->nextEdge = NULL; // now the last edge


        // For everyguard inside the polygon
        for(k=0; k<insideGuardArray[0].isNull; k++){
            printf("guards: %d \n", insideGuardArray[0].isNull);
            Vertice curGuard = insideGuardArray[k];

            // Replicate polygon array.
            Vertice* polygonReplica = calloc(polygonArray[0].isNull, sizeof(Vertice));
            copy_vertice_array(polygonReplica, polygonArray);

            // Update polar angle of all polygon vertex to guard k
            for (int i = 0; i < polygonReplica[0].isNull; i++){
                polygonReplica[i].angle_to_guard = polar_angle(polygonReplica[i].x, polygonReplica[i].y, curGuard.x, curGuard.y)+M_PI;                
            }
            // Sort by polar angle to guard k
            qsort(polygonReplica, polygonReplica[0].isNull, sizeof(struct vertice), &md_comparator);

            // Construct a new visibility Polygon
            Vertice* visiblePolygon = calloc( 1000, sizeof(struct vertice));
            int counter = 0;

            // Finding visible point
            for(int i=0; i<polygonReplica[0].isNull; i++){   // For every vertices in sorted array
                Vertice curVert = polygonReplica[i];

                double distance = sqrt(pow(curGuard.x - curVert.x , 2) + pow(curGuard.y - curVert.y , 2)); // distance from guard to vertex


                // Cast a ray from guard to vertex
                // as a line segment from guard to temp
                // y = mx + c
                double m = (curVert.y - curGuard.y)/(curVert.x - curGuard.x);
                double c = curGuard.y - m*curGuard.x;
                Vertice* temp = calloc(1, sizeof(struct vertice));
                temp->x =1000; // very far away
                temp->y = m*(temp->x) + c;

                // Get first intersection point and on which edge
                
                Edge* min_edge = first;
                double minX = curVert.x, minY = curVert.y;
                Edge* curEdge = first;
                double min_distance = 9999;
                while(curEdge->nextEdge != NULL){ // For every edge in original polygon
                    double x0, y0;
                    int isEndPoint = get_line_intersection(curGuard.x, curGuard.y,temp->x, temp->y, curEdge->strt->x, curEdge->strt->y, curEdge->fins->x, curEdge->fins->y, &x0, &y0);
                    if(isEndPoint == -1) { // If x0 y0 belongs to either the guard  ( guard sit on polygon vertice/edge) or the Edge skip 
                        curEdge= curEdge->nextEdge;
                        continue;
                    }
                    double temp_distance = sqrt(pow(curGuard.x - x0, 2 ) + pow(curGuard.y - y0, 2 )); //dist from intersect point to guard
                    if( fabs(temp_distance-min_distance)<EPSI ){
                        min_edge = curEdge;
                        min_distance = temp_distance;
                        minX = x0, minY = y0;
                    }
                    curEdge= curEdge->nextEdge;
                }
                // after while min_edge = nearest wall hit by ray
                // minX minY now the coordinate of nearest intersecting
                // if isEndPoint == -1 for all edges, minX minY is the current vertex
                        

                // consider min_distance = distance
                    // can't happen due to isEndPoint = -1 thus skip

                if ( min_distance < distance - EPSI ){ 
                    // collision inside line segment
                    
                    Vertice* new = malloc(sizeof(struct vertice));
                    new->x = minX, new->y = minY;
                    

                    // Split the edge
                    Edge* new_edge = malloc(sizeof(struct edge));
                    new_edge->strt = new;
                    new_edge->fins = min_edge->fins;
                    new_edge->nextEdge = min_edge->nextEdge;
                    new_edge->prevEdge = min_edge;
                    min_edge->nextEdge = new_edge;
                    min_edge->fins = new;

                    // add the collision vertex
                    visiblePolygon[counter] = *new;

                }
                if( min_distance > distance + EPSI ){ 
                    //If Outside line segment (not on the end points of the edge)
                    // Divide the intersect edge to 2 new edge
                    
                    // add the collision vertex
                    Vertice* new = malloc(sizeof(struct vertice));
                    new->x = minX, new->y = minY;
                    

                    // Split the edge
                    Edge* new_edge = malloc(sizeof(struct edge));
                    new_edge->strt = new;
                    new_edge->fins = min_edge->fins;
                    new_edge->nextEdge = min_edge->nextEdge;
                    new_edge->prevEdge = min_edge;
                    min_edge->nextEdge = new_edge;
                    min_edge->fins = new;

                    // add the collision vertex
                    visiblePolygon[counter] = *new;
                    counter++; // in this case counter is increment twice
                    // then add the curVertex
                    visiblePolygon[counter] = curVert;

                }
                if (fabs(minX - curVert.x) < EPSI && fabs(minY - curVert.y)< EPSI){ // no edge collision
                    // add the curVert
                    visiblePolygon[counter] = curVert;    
                }
                counter++;
                // Add vertex to visibility polygon
                
                free(temp);
            }
            free(polygonReplica);
            free(visiblePolygon);
        }
    } 
}
