
#!/bin/bash 

index=0 
> makefile 
for myfile in AUX/*.c

do 
    string=${myfile%.*}
    if [ $index -eq 0 ] 
    then 
        printf "${string##*/} " >> makefile 
    else 
        ((index++)) 
        printf "${string##*/} " >> makefile 
    fi 
done

echo '#' >> makefile
for myfile in *.c
do 
    string=${myfile%.c}
    if [ $string = "cs" ]
    then
        continue
    else
        if [ $index -eq 0 ] 
        then 
            printf "${string##*/} " >> makefile 
        else 
            ((index++)) 
            printf "${string##*/} " >> makefile 
        fi 
    fi
done
echo '#' >> makefile



