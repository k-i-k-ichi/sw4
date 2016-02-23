#include <stdio.h>

int charcount( FILE *const fin )
{
    int c, count;
    count = 0;
    for( ;; )
    {
        c = fgetc( fin );
        if( c == EOF || c == '\n' )
            break;
        ++count;
    }

    return count;
}

int main(){
    FILE* fp;
    fp = fopen("check.pol");
    for(int main=0; main<20; main++){
        int lineCharCount = charcount(fp);
            
    } 
}
