/**
 * @file proj2.c
 * 
 * @author Ela Fedorová (xfedor17)
 * @brief Druhý projekt do predmetu IOS - synchronizácia
 * @date 2022-04-29
 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <semaphore.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/sem.h>
#include <sys/mman.h>
#include <sys/wait.h>

//premenne pre zdielanu pamat
int *serial_number = NULL;
int *oxygens_to_bond = NULL;
int *hydrogens_to_bond = NULL;
int *molecule_id = NULL;
int *max_molecules = NULL;       //maximalny pocet molekul, ktore mozu byt vytvorene
bool *all_molecules = NULL;      //na urcenie, ci uz su vytvorene vsetky molekuly, ktore mozu byt vytvorene

//semafory
sem_t *oxygen_queue = NULL;
sem_t *hydrogen_queue = NULL;
sem_t *mutex = NULL;
sem_t *write_msg = NULL;
sem_t *bond = NULL;
sem_t *oxy_created = NULL;         //proces kyslik vytvoril molekulu - sprava pre vodik
sem_t *hydro_creating = NULL;      //proces vodik vytvara molekulu - sprava pre kyslik
sem_t *hydro_created = NULL;       //proces vodik vytvoril molekulu - sprava pre kyslik

FILE *file_out;

/*
* Funkcia initialize() na otvorenie vystupneho suboru, inicializaciu zdielanej pamate a semaforov
*/
void initialize()
{
    if((file_out = fopen("proj2.out", "w")) == NULL)
    {
        fprintf(stderr, "Error: file opening failed.\n");
        exit(EXIT_FAILURE);
    }

    //inicializacia zdielanej pamate
    if((serial_number = mmap(NULL, sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0)) == MAP_FAILED ||
    (oxygens_to_bond = mmap(NULL, sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0)) == MAP_FAILED ||
    (hydrogens_to_bond = mmap(NULL, sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0)) == MAP_FAILED ||
    (molecule_id = mmap(NULL, sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0)) == MAP_FAILED ||
    (max_molecules = mmap(NULL, sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0)) == MAP_FAILED ||
    (all_molecules = mmap(NULL, sizeof(bool), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0)) == MAP_FAILED)
    {
        fprintf(stderr, "Error: shared memory initialization failed.\n");
        exit(EXIT_FAILURE);
    }

    //inicializacia semaforov
    if((oxygen_queue = mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0)) == MAP_FAILED ||
    (hydrogen_queue = mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0)) == MAP_FAILED ||
    (mutex = mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0)) == MAP_FAILED ||
    (write_msg = mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0)) == MAP_FAILED ||
    (bond = mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0)) == MAP_FAILED ||
    (oxy_created = mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0)) == MAP_FAILED ||
    (hydro_creating = mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0)) == MAP_FAILED ||
    (hydro_created = mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0)) == MAP_FAILED)
    {
        fprintf(stderr, "Error: shared memory for semaphores - initialization failed.\n");
        exit(EXIT_FAILURE);
    }

    if((sem_init(oxygen_queue, 1, 0)) == -1 || (sem_init(hydrogen_queue, 1, 0)) == -1 || (sem_init(mutex, 1, 1)) == -1 || (sem_init(write_msg, 1, 1)) == -1 ||
    (sem_init(bond, 1, 1)) == -1 || (sem_init(oxy_created, 1, 0)) == -1 || (sem_init(hydro_creating, 1, 0)) == -1 || (sem_init(hydro_created, 1, 0)) == -1)
    {
        fprintf(stderr, "Error: semaphore initialization failed\n");
        exit(EXIT_FAILURE);
    }
}

/*
* Funkcia clean() na zatvorenie vystupneho suboru, odmapovanie zdielanej pamate a znicenie semaforov
*/
void clean()
{
    fclose(file_out);

    //odmapovanie zdielanej pamate
    if((munmap((serial_number), sizeof(serial_number))) == -1 ||
    (munmap((oxygens_to_bond), sizeof(oxygens_to_bond))) == -1 ||
    (munmap((hydrogens_to_bond), sizeof(hydrogens_to_bond))) == -1 ||
    (munmap((molecule_id), sizeof(molecule_id))) == -1 ||
    (munmap((max_molecules), sizeof(max_molecules))) == -1 ||
    (munmap((all_molecules), sizeof(all_molecules))) == -1)
    {
        fprintf(stderr, "Error: shared memory deinitialization failed.\n");
        exit(EXIT_FAILURE);
    }

    //ukoncenie semaforov
    if((sem_destroy(oxygen_queue)) == -1 || (sem_destroy(hydrogen_queue)) == -1 || (sem_destroy(mutex)) == -1 || (sem_destroy(write_msg)) == -1 ||
    (sem_destroy(bond)) == -1 || (sem_destroy(oxy_created)) == -1 || (sem_destroy(hydro_creating)) == -1 || (sem_destroy(hydro_created)) == -1)
    {
        fprintf(stderr, "Error: Semaphore destruction failed.\n");
        exit(EXIT_FAILURE);
    }

    exit(EXIT_SUCCESS);
}

/*
* Funkcia oxygen_process() na cely proces kyslik
* Parametre: id - id kyslika (podla poradia procesov), TI - maximalny cas cakania, kym sa kyslik zaradi do fronty,
*            TB - maximalny cas nutny pre vytvorenie jednej molekuly
*/
void oxygen_process(int id, int TI, int TB) 
{   
    srand(getpid()^time(NULL));
    if((*max_molecules) == 0) {             //ak nie je mozne vytvorit ziadne molekuly, rovno sa otvori fronta kyslikov
        sem_post(oxygen_queue);
    }
    id++;
    sem_wait(write_msg);
        fprintf(file_out, "%d: O %d: started\n", ++(*serial_number), id);
        fflush(file_out);
    sem_post(write_msg);

    usleep((rand() % (TI + 1)) * 1000);

    sem_wait(write_msg);
        fprintf(file_out, "%d: O %d: going to queue\n", ++(*serial_number), id);
        fflush(file_out);
    sem_post(write_msg);

    sem_wait(bond);
    sem_wait(mutex);
        (*oxygens_to_bond)++;
        if((*hydrogens_to_bond) >= 2)       //dostatok cakajucich vodikov, moze sa vytvarat molekula
        {
            (*molecule_id)++;

            sem_post(hydrogen_queue);       //dvakrat signal na uvolnenie cakajucich vodikov    
            sem_post(hydrogen_queue);
            (*hydrogens_to_bond) -= 2;      //odstranenie dvoch vodikov z fronty cakajucich

            sem_post(oxygen_queue);         //signal na uvolnenie cakajuceho kyslika
            (*oxygens_to_bond)--;           //odstranenie kyslika z fronty cakajucich

            sem_post(mutex);
        }
        else                                //nie je dostatok vodikov, nebude sa vytvarat molekula
        {
            sem_post(bond);
            sem_post(mutex);
        }

        sem_wait(oxygen_queue);                             //cakanie na signal o uvolneni kyslika na vytvaranie molekuly
        if((*max_molecules) == 0 || (*all_molecules))       //ak nie je mozne vytvorit molekulu
        {
            sem_wait(write_msg);
                fprintf(file_out, "%d: O %d: not enough H\n", ++(*serial_number), id);
                fflush(file_out);
            sem_post(write_msg);
            sem_post(oxygen_queue);          //uvolnenie zvysnych kyslikov a ukoncenie procesu
            exit(EXIT_SUCCESS);
        }

        sem_wait(write_msg);
            fprintf(file_out, "%d: O %d: creating molecule %d\n", ++(*serial_number), id, (*molecule_id));
            fflush(file_out);
        sem_post(write_msg);

        usleep((rand() % (TB + 1)) * 1000);         //simuluje dobu vytvarania molekuly

        sem_wait(hydro_creating);                   //caka na signal od dvoch vodikov, ze vytvaraju molekulu - skoncili vypis creating
        sem_wait(hydro_creating);
        sem_post(oxy_created);                      //dava signal dvom vodikom, ze molekula bola vytvorena
        sem_post(oxy_created);
        
        sem_wait(write_msg);
            fprintf(file_out, "%d: O %d: molecule %d created\n", ++(*serial_number), id, (*molecule_id));
            fflush(file_out);
        sem_post(write_msg);

        sem_wait(hydro_created);                    //caka na signal od dvoch vodikov, ze skoncili vypis created
        sem_wait(hydro_created);

        sem_wait(mutex);
            if((*molecule_id) == (*max_molecules))      //ak uz su vytvorene vsetky molekuly
            {
                (*all_molecules) = true;
                sem_post(oxygen_queue);                 //uvolnenie zvysnych kyslikov a vodikov
                sem_post(hydrogen_queue);
            }
        sem_post(mutex);

    sem_post(bond);
    exit(EXIT_SUCCESS);
}

/*
* Funkcia hydrogen_process() na cely proces vodik
* Parametre: id - id vodika (podla poradia procesov), TI - maximalny cas cakania, kym sa vodik zaradi do fronty
*/
void hydrogen_process(int id, int TI)
{   
    srand(getpid()^time(NULL));
    if((*max_molecules) == 0) {                     //ak nie je mozne vytvorit ziadne molekuly, rovno sa otvori fronta vodikov
        sem_post(hydrogen_queue);
    }

    id++;
    sem_wait(write_msg);
        fprintf(file_out, "%d: H %d: started\n", ++(*serial_number), id);
        fflush(file_out);
    sem_post(write_msg);

    usleep((rand() % (TI + 1)) * 1000);
    
    sem_wait(write_msg);
        fprintf(file_out, "%d: H %d: going to queue\n", ++(*serial_number), id);
        fflush(file_out);
    sem_post(write_msg);

    sem_wait(bond);
    sem_wait(mutex);
        (*hydrogens_to_bond)++;
        if((*hydrogens_to_bond) >= 2 && (*oxygens_to_bond) >= 1)        //dostatok cakajucich vodikov a kyslikov, moze sa vytvorit molekula
        {
            (*molecule_id)++;

            sem_post(hydrogen_queue);               //dvakrat signal na uvolnenie cakajucich vodikov
            sem_post(hydrogen_queue);
            (*hydrogens_to_bond) -= 2;              //odstranenie dvoch vodikov z fronty cakajucich

            sem_post(oxygen_queue);                 //signal na uvolnenie cakajuceho kyslika
            (*oxygens_to_bond)--;                   //odstranenie kyslika z fronty cakajucich

            sem_post(mutex);
        }
        else                                        //nie je dostatok vodikov alebo kyslikov, nebude sa vytvarat molekula
        {
            sem_post(bond);
            sem_post(mutex);
        }

        sem_wait(hydrogen_queue);                           //cakanie na signal o uvolneni prveho vodika na vytvaranie molekuly
        if((*max_molecules) == 0 || (*all_molecules))       //ak nie je mozne vytvorit molekulu
        {
            sem_wait(write_msg);
                fprintf(file_out, "%d: H %d: not enough O or H\n", ++(*serial_number), id);
                fflush(file_out);
            sem_post(write_msg);
            sem_post(hydrogen_queue);          //uvolnenie zvysnych vodikov a ukoncenie procesu
            exit(EXIT_SUCCESS);
        }

        sem_wait(write_msg);
            fprintf(file_out, "%d: H %d: creating molecule %d\n", ++(*serial_number), id, (*molecule_id));
            fflush(file_out);
        sem_post(write_msg);
        
        sem_post(hydro_creating);              //posiela signal kysliku, ze zacal vytvarat molekulu - vypisal creating
        sem_wait(oxy_created);                 //caka na signal od kyslika, ze je molekula dokoncena

        sem_wait(write_msg);
            fprintf(file_out, "%d: H %d: molecule %d created\n", ++(*serial_number), id, (*molecule_id));
            fflush(file_out);
            sem_post(hydro_created);           //posiela signal kysliku, ze ukoncil vypis created
        sem_post(write_msg);
    
    exit(EXIT_SUCCESS);
}


/*
* Funkcia arg_check() kontroluje spravnost zadanych argumentov
* Parametre: argc - pocet argumentov, argv - ukazatel na pole argumentov
*/
void arg_check(int argc, char **argv)
{
    if(argc != 5)
    {
        fprintf(stderr, "Incorrect number of arguments! Run the programm with: ./proj2 NO NH TI TB\n");
        exit(EXIT_FAILURE);
    }

    for(int i=1; i<argc; i++)
    {
        char *arg = argv[i];
        for(int j=0; arg[j] != '\0'; j++)
        {
            if(!isdigit(arg[j]))
            {
                fprintf(stderr, "Invalid argument! All arguments have to be positive integers or 0.\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    if(atoi(argv[3]) > 1000 || atoi(argv[4]) > 1000)
    {
        fprintf(stderr, "Invalid value! Values of parameters TI and TB must be in range <0;1000>.\n");
        exit(EXIT_FAILURE);
    }

    if(atoi(argv[1]) == 0 || atoi(argv[2]) == 0)
    {
        fprintf(stderr, "Invalid value! Values of parameters NO and NH must be greater than 0.\n");
        exit(EXIT_FAILURE);
    }
}

/*
* Funkcia main()
*/
int main(int argc, char **argv)
{
    arg_check(argc, argv);

    int NO = atoi(argv[1]);
    int NH = atoi(argv[2]);
    int TI = atoi(argv[3]);
    int TB = atoi(argv[4]);

    initialize();

    *serial_number = 0;
    *oxygens_to_bond = 0;
    *hydrogens_to_bond = 0;
    *molecule_id = 0;
    *all_molecules = false;

    //maximalny mozny pocet molekul na vytvorenie
    if((NH / 2) <= NO) { *max_molecules = (NH /2); }
    else { *max_molecules = NO; }

    //hydrogen process generator
    for(int i = 0; i < NH; i++)
    {
        pid_t pid = fork();
        if(pid == 0)
        {
            hydrogen_process(i, TI);
            exit(EXIT_SUCCESS);
        }
        else if(pid < 0)
        {
            fprintf(stderr, "Error: fork failed.\n");
            clean();
            exit(EXIT_FAILURE);       
        }
    }  
    
    //oxygen process generator
    for(int i = 0; i < NO; i++)
    {
        pid_t pid = fork();
        if(pid == 0)
        {
            oxygen_process(i, TI, TB);
            exit(EXIT_SUCCESS);
        }
        else if(pid < 0)
        {
            fprintf(stderr, "Error: fork failed.\n");                
            clean();
            exit(EXIT_FAILURE);       
        }
    }
    
    while(wait(NULL) > 0);              //cakanie, kym skoncia child procesy

    clean();
    exit(EXIT_SUCCESS);
}
