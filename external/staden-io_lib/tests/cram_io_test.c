#define CRAM_IO_TEST
#include <io_lib/scram.h>
#include <assert.h>
#include <string.h>

int main(int argc, char *argv[])
{
    int i = 0;
    for ( i = 1; i < argc; ++i ) {
        FILE * fp = fopen(argv[i],"rb");
        cram_fd * cramfd = NULL;
        char * Ba = NULL;
        char * Bb = NULL;
        size_t la = 0;
        size_t lb = 0;
        size_t o = 0;
        size_t p = 0;
        int r = -1;
        char linebuf0[32];
        char linebuf1[32];

        if ( ! fp ) {
            fprintf(stderr,"Cannot open file %s\n",argv[i]);
            goto cleanup;
        }
        
        cramfd = cram_io_open(argv[i],"rc","rb");
        if ( ! cramfd )
            goto cleanup;

        /* compare file sizes by seeking to end of file */
        r = fseek(fp,0,SEEK_END);
        assert ( r == 0 );
        la = ftello(fp);
        
        r = CRAM_IO_SEEK(cramfd,0,SEEK_END);
        assert ( r == 0 );
        lb = CRAM_IO_TELLO(cramfd);
        
        assert ( la == lb );
        
        Ba = (char *)malloc(la);
        if ( ! Ba )
	    goto cleanup;
        Bb = (char *)malloc(lb);
        if ( ! Bb )
	    goto cleanup;
        
        /* seek to start position and read file via getc type calls */
        for ( o = 0; o <= la; ++o ) {
            r = fseeko(fp,o,SEEK_SET);
            assert ( r == 0 );
            r = CRAM_IO_SEEK(cramfd,o,SEEK_SET);
            assert ( r == 0 );

            if ( o % 1024 == 0 ) {
                fprintf(stderr,"%s/%d/%d\n",argv[i],(int)o, (int)la);
	    }
            
            for ( p = o; p < la; ++p ) {
                int const c1 = getc(fp);
                int const c2 = CRAM_IO_GETC(cramfd);
                assert ( c1 == c2 );
            }
            
            assert ( getc(fp) == EOF );
            assert ( CRAM_IO_GETC(cramfd) == EOF );
        }

        /* seek to start position, read up to 16 bytes via getc and the rest via fread() */
        for ( o = 0; o <= la; ++o ) {
            size_t rest = 0;
            int r0 = -1;
            int r1 = -1;
            
            r = fseeko(fp,o,SEEK_SET);
            assert ( r == 0 );
            r = CRAM_IO_SEEK(cramfd,o,SEEK_SET);
            assert ( r == 0 );

            if ( o % 1024 == 0 ) {
                fprintf(stderr,"%s/%d/%d\n",argv[i],(int)o, (int)la);
	    }
	    
	    for ( p = o; p < la && p < o+16; ++p ) {
                int const c1 = getc(fp);
                int const c2 = CRAM_IO_GETC(cramfd);
                assert ( c1 == c2 );	    
	    }
	    
	    assert ( p <= la );
	    rest = la-p;
	    	    
	    r0 = fread(Ba,rest,1,fp);
	    assert ( rest == 0 || r0 == 1 );
	    r1 = CRAM_IO_READ(Bb,rest,1,cramfd);
	    assert ( rest == 0 || r1 == 1 );
	    
	    assert ( r0 == r1 );
	    
	    for ( p = 0; p < rest; ++p ) {
	        assert ( Ba[p] == Bb[p] );
	    }

            assert ( getc(fp) == EOF );
            assert ( CRAM_IO_GETC(cramfd) == EOF );
	}

        /* run fgets type test */
        r = fseeko(fp,0,SEEK_SET);
        assert ( r == 0 );
        r = CRAM_IO_SEEK(cramfd,0,SEEK_SET);
        assert ( r == 0 );
        while ( fgets(&linebuf0[0],sizeof(linebuf0),fp) ) {
            assert ( CRAM_IO_FGETS(&linebuf1[0],sizeof(linebuf1),cramfd) == &linebuf1[0] );
            assert ( strlen(&linebuf0[0]) == strlen(&linebuf1[0]) );
            assert ( strcmp(&linebuf0[0],&linebuf1[0]) == 0 );
        }
        assert ( CRAM_IO_FGETS(&linebuf1[0],sizeof(linebuf1),cramfd) == NULL );
        
        cleanup:
        if ( Bb ) {
            free(Bb);
            Bb = NULL;
        }
        if ( Ba ) {
            free(Ba);
            Ba = NULL;
        }
        if ( cramfd ) {
            cram_io_close(cramfd,NULL);
            cramfd = NULL;
        }
        if ( fp ) {
            fclose(fp);
            fp = NULL;
        }
    }

    return 0;
}
