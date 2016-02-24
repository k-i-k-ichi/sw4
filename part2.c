#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct vertice{
    double x;
    double y;
}Vertice;


Vertice* parse(char* buffer){
    int memberCount=0;
    int i;
    while(buffer[i]!=NULL){
        if(buffer[i] == '(') memberCount++;
        i++;
    }
    Vertice* verticeArray = malloc(memberCount*sizeof(struct vertice));
    
    int pointer=1;
    
    for(int i=0; i< memberCount;i++){
        char doubleBuffer[100];
        int initialPointer = pointer;
        
        while(buffer[pointer]!= ','){ // x coordinate
            doubleBuffer[pointer-initialPointer] = buffer[pointer];
            pointer++;
        }
        pointer+=2; // move past char ','
        
        
        char doubleBuffer2[100];
        initialPointer = pointer;
        while(buffer[pointer]!=')'){ // y coordinate
            doubleBuffer2[pointer-initialPointer] = buffer[pointer];
            pointer++;
        }
        
        while(buffer[pointer]!='(') pointer++; // seek to next '('
        
        pointer++; // skip ( to next number)
        
        Vertice* new = malloc(sizeof(struct vertice));
        new->x = atof(doubleBuffer);
        new->y = atof(doubleBuffer2);
        
        verticeArray[i] = *new;
    }
    
    return verticeArray;
}


int main(){
    FILE* fp;
    fp = fopen("check.pol", "r");
    for(int main=0; main<1; main++){
        char* buffer = malloc(25000*sizeof(char));
        fgets(buffer, 25000, fp);
        
        puts(buffer);
        int i =3;
        
        char* polygon = malloc(20000*sizeof(char));
        while(buffer[i] != ';'){
             polygon[i] = buffer[i];
             i++;
        }
        
        polygon[i+1]='\0';
        
        puts(polygon);
        
        i+= 2;
        char* guard = malloc(8000*sizeof(char));
        while(buffer[i] != '\n'){
                guard[i] = buffer[i];
                i++;
        }
        guard[i+1]='\0';
        free(buffer);
         
    } 
}
