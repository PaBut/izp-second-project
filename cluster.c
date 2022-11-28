/**Pavlo Butenko
 * xbuten00
 * Kostra programu pro 2. projekt IZP 2022/23
 * 
 * Jednoducha shlukova analyza: 2D nejblizsi soused.
 * Single linkage
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h> // sqrtf
#include <limits.h> // INT_MAX
#include <string.h>
#include <stdbool.h>
//#include <errno.h>

/*****************************************************************
 * Ladici makra. Vypnout jejich efekt lze definici makra
 * NDEBUG, napr.:
 *   a) pri prekladu argumentem prekladaci -DNDEBUG
 *   b) v souboru (na radek pred #include <assert.h>
 *      #define NDEBUG
 */
#ifdef NDEBUG
#define debug(s)
#define dfmt(s, ...)
#define dint(i)
#define dfloat(f)
#else

// vypise ladici retezec
#define debug(s) printf("- %s\n", s)

// vypise formatovany ladici vystup - pouziti podobne jako printf
#define dfmt(s, ...) printf(" - "__FILE__":%u: "s"\n",__LINE__,__VA_ARGS__)

// vypise ladici informaci o promenne - pouziti dint(identifikator_promenne)
#define dint(i) printf(" - " __FILE__ ":%u: " #i " = %d\n", __LINE__, i)

// vypise ladici informaci o promenne typu float - pouziti
// dfloat(identifikator_promenne)
#define dfloat(f) printf(" - " __FILE__ ":%u: " #f " = %g\n", __LINE__, f)

#endif

/*****************************************************************
 * Deklarace potrebnych datovych typu:
 *
 * TYTO DEKLARACE NEMENTE
 *
 *   struct obj_t - struktura objektu: identifikator a souradnice
 *   struct cluster_t - shluk objektu:
 *      pocet objektu ve shluku,
 *      kapacita shluku (pocet objektu, pro ktere je rezervovano
 *          misto v poli),
 *      ukazatel na pole shluku.
 */

struct obj_t {
    int id;
    float x;
    float y;
};

struct cluster_t {
    int size;
    int capacity;
    struct obj_t *obj;
};

/*****************************************************************
 * Deklarace potrebnych funkci.
 *
 * PROTOTYPY FUNKCI NEMENTE
 *
 * IMPLEMENTUJTE POUZE FUNKCE NA MISTECH OZNACENYCH 'TODO'
 *
 */

/*
 Inicializace shluku 'c'. Alokuje pamet pro cap objektu (kapacitu).
 Ukazatel NULL u pole objektu znamena kapacitu 0.
*/
void init_cluster(struct cluster_t *c, int cap)
{
    assert(c != NULL);
    assert(cap >= 0);

    // TODO
    if((c->obj = malloc(cap * sizeof(struct obj_t))) == NULL){
        return;
    }
    c->capacity = cap;
    c->size = 0;
}

/*
 Odstraneni vsech objektu shluku a inicializace na prazdny shluk.
 */
void clear_cluster(struct cluster_t *c)
{
    // TODO
    free(c->obj);
    c->size = 0;
    c->capacity = 0;
}

/// Chunk of cluster objects. Value recommended for reallocation.
const int CLUSTER_CHUNK = 10;

/*
 Zmena kapacity shluku 'c' na kapacitu 'new_cap'.
 */
struct cluster_t *resize_cluster(struct cluster_t *c, int new_cap)
{
    // TUTO FUNKCI NEMENTE
    assert(c);
    assert(c->capacity >= 0);
    assert(new_cap >= 0);

    if (c->capacity >= new_cap)
        return c;

    size_t size = sizeof(struct obj_t) * new_cap;

    void *arr = realloc(c->obj, size);
    if (arr == NULL)
        return NULL;

    c->obj = (struct obj_t*)arr;
    c->capacity = new_cap;
    return c;
}

/*
 Prida objekt 'obj' na konec shluku 'c'. Rozsiri shluk, pokud se do nej objekt
 nevejde.
 */
void append_cluster(struct cluster_t *c, struct obj_t obj)
{
    // TODO
    if(c->size + 1 > c->capacity){
        if(resize_cluster(c, c->size+1) == NULL){
            return;
        }
    }
    c->obj[c->size++] = obj;
}

/*
 Seradi objekty ve shluku 'c' vzestupne podle jejich identifikacniho cisla.
 */
void sort_cluster(struct cluster_t *c);

/*
 Do shluku 'c1' prida objekty 'c2'. Shluk 'c1' bude v pripade nutnosti rozsiren.
 Objekty ve shluku 'c1' budou serazeny vzestupne podle identifikacniho cisla.
 Shluk 'c2' bude nezmenen.
 */
void merge_clusters(struct cluster_t *c1, struct cluster_t *c2)
{
    assert(c1 != NULL);
    assert(c2 != NULL);

    // TODO
    //Checking here if merging can fit in allocated memory for objects in c1
    if(c1->size + c2->size > c1->capacity){
        if(resize_cluster(c1, c1->size + c2->size) == NULL){
            return;
        }
    }
    for(int i = 0; i < c2->size; i++){
        c1->obj[c1->size + i] = c2->obj[i];
    }
    c1->size+=c2->size;
    sort_cluster(c1);
}

/**********************************************************************/
/* Prace s polem shluku */

/*
 Odstrani shluk z pole shluku 'carr'. Pole shluku obsahuje 'narr' polozek
 (shluku). Shluk pro odstraneni se nachazi na indexu 'idx'. Funkce vraci novy
 pocet shluku v poli.
*/
int remove_cluster(struct cluster_t *carr, int narr, int idx)
{
    assert(idx < narr);
    assert(narr > 0);

    // TODO
    clear_cluster(&carr[idx]);

    for(int i = idx; i < narr - 1; i++){
        carr[i] = carr[i+1];
    }
    return narr-1;
}

/*
 Pocita Euklidovskou vzdalenost mezi dvema objekty.
 */
float obj_distance(struct obj_t *o1, struct obj_t *o2)
{
    assert(o1 != NULL);
    assert(o2 != NULL);

    // TODO
    return sqrt(pow(o1->x - o2->x, 2) + pow(o1->y - o2->y, 2));
}

/*
 Pocita vzdalenost dvou shluku.
*/
float cluster_distance(struct cluster_t *c1, struct cluster_t *c2)
{
    assert(c1 != NULL);
    assert(c1->size > 0);
    assert(c2 != NULL);
    assert(c2->size > 0);

    // TODO
    float min_distance = obj_distance(&c1->obj[0], &c2->obj[0]);
    for(int i = 0; i < c1->size; i++){
        for(int j = 0; j < c2->size; j++){
            float tmp_distance = obj_distance(&c1->obj[i], &c2->obj[j]);
            if(tmp_distance < min_distance){
                min_distance = tmp_distance;
            }       
        }
    }
    return min_distance;
}

/*
 Funkce najde dva nejblizsi shluky. V poli shluku 'carr' o velikosti 'narr'
 hleda dva nejblizsi shluky. Nalezene shluky identifikuje jejich indexy v poli
 'carr'. Funkce nalezene shluky (indexy do pole 'carr') uklada do pameti na
 adresu 'c1' resp. 'c2'.
*/
void find_neighbours(struct cluster_t *carr, int narr, int *c1, int *c2)
{
    assert(narr > 0);
    // TODO

    float min_distance = cluster_distance(&carr[0], &carr[1]);
    for(int i = 0; i < narr; i++){
        for(int j = i+1; j < narr; j++){
            int tmp_distance = cluster_distance(&carr[i], &carr[j]);
            if(tmp_distance <= min_distance){
                min_distance = tmp_distance;
                *c1 = i;
                *c2 = j;
            }
        }
    }
}

// pomocna funkce pro razeni shluku
static int obj_sort_compar(const void *a, const void *b)
{
    // TUTO FUNKCI NEMENTE
    const struct obj_t *o1 = (const struct obj_t *)a;
    const struct obj_t *o2 = (const struct obj_t *)b;
    if (o1->id < o2->id) return -1;
    if (o1->id > o2->id) return 1;
    return 0;
}

/*
 Razeni objektu ve shluku vzestupne podle jejich identifikatoru.
*/
void sort_cluster(struct cluster_t *c)
{
    // TUTO FUNKCI NEMENTE
    qsort(c->obj, c->size, sizeof(struct obj_t), &obj_sort_compar);
}

/*
 Tisk shluku 'c' na stdout.
*/
void print_cluster(struct cluster_t *c)
{
    // TUTO FUNKCI NEMENTE
    for (int i = 0; i < c->size; i++)
    {
        if (i) putchar(' ');
        printf("%d[%g,%g]", c->obj[i].id, c->obj[i].x, c->obj[i].y);
    }
    putchar('\n');
}

void free_clusters(struct cluster_t** arr, int len);
int return_error(char* msg, int return_value);

/*
 Ze souboru 'filename' nacte objekty. Pro kazdy objekt vytvori shluk a ulozi
 jej do pole shluku. Alokuje prostor pro pole vsech shluku a ukazatel na prvni
 polozku pole (ukalazatel na prvni shluk v alokovanem poli) ulozi do pameti,
 kam se odkazuje parametr 'arr'. Funkce vraci pocet nactenych objektu (shluku).
 V pripade nejake chyby uklada do pameti, kam se odkazuje 'arr', hodnotu NULL.
*/
int load_clusters(char *filename, struct cluster_t **arr)
{

    //assert(arr != NULL);

    // TODO
    FILE* file;
    file = fopen(filename, "r");
    if(file == NULL){
        fprintf(stderr, "Argument error\n%s doesn\'t exist\n", filename);
        return 0;
    }
    char count_txt[31];

    fgets(count_txt, 30, file);
    //Checking here if a file in a proper format(must contain count= in the beginning)
    if(strcmp(strtok(count_txt, "="), "count") != 0){
        fclose(file);
        return return_error("Argument error\nArgument file isn\'t in a proper format", 0);
    }
    int count = atoi(strtok(NULL, "="));
    //Checking here if count argument in file is an integer higher than 0
    if(count <= 0){
        fclose(file);
        return return_error("Argument file error\nThe value of [count] must be integer higher than 0", 0);
    }

    *arr = malloc(count * sizeof(struct cluster_t));
    if(arr == NULL){
        fclose(file);
        return return_error("Allocation error", 0);
    }
    int i = 0;
    //Setting file values to clusters
    for(; i < count; i++){
        struct obj_t obj;
        int scan_res = fscanf(file, "%d %f %f", &obj.id, &obj.x, &obj.y);
        if(scan_res == EOF){
            fclose(file);
            fprintf(stderr, "Argument file error\nCount of objects given \
in the beginning of the file is higher than initial count of objects in the file");
            return i;
        }
        if(scan_res != 3){
            fclose(file);
            free_clusters(arr, i);
            fprintf(stderr, "Argument file error\nProperies of objects isn\'t proper, they must be int type");
            return 0;
        }
        if(obj.x < 0 || obj.x >1000 || obj.y < 0 || obj.y > 1000){
            fclose(file);
            free_clusters(arr, i);
            return return_error("Argument file error\nThe value of coordinates \
of some cluster must be in range between 0 and 1000 included", 0);
        }
        init_cluster(&(*arr)[i], 1);
        append_cluster(&(*arr)[i], obj);
        (*arr)[i].size = 1;
        /*init_cluster(&(*arr)[i], 1);
        fscanf(file, "%d %f %f", &(*arr)[i].obj->id, &(*arr)[i].obj->x, &(*arr)[i].obj->y);
        (*arr)[i].size = 1;*/
    }

    fclose(file);
    //Checking if all objects have unique id
    for(int j = 0; j < i; j++){
        for(int k = j + 1; k < i; k++){
            if((*arr)[j].obj[0].id == (*arr)[k].obj[0].id){
                free_clusters(arr, i);
                fprintf(stderr, "Argument file error\nID in a cluster isn\'t unique");
                return 0;
            }
        }
    }
    return i;
}

/*
 Tisk pole shluku. Parametr 'carr' je ukazatel na prvni polozku (shluk).
 Tiskne se prvnich 'narr' shluku.
*/
void print_clusters(struct cluster_t *carr, int narr)
{
    printf("Clusters:\n");
    for (int i = 0; i < narr; i++)
    {
        printf("cluster %d: ", i);
        print_cluster(&carr[i]);
    }
}


//Function to deallocate array of clusters
void free_clusters(struct cluster_t** arr, int len){
    for(int i = 0; i < len; i++){
        clear_cluster(&(*arr)[i]);
    }
    free(*arr);
}

//Function to return error with message to stderr
int return_error(char* msg, int return_value){
    fprintf(stderr, "%s\n", msg);
    return return_value;
}

int main(int argc, char *argv[])
{
    struct cluster_t *clusters;

    //Checking if arguments given to program are proper
    if(!(argc == 2 || argc == 3)){
        return return_error("Argument error\nWrong amount of arguments", EXIT_FAILURE);
    }

    int count = 1;

    if(argc == 3){
        int remainder = 0;
        if(sscanf(argv[2], "%d.%d", &count, &remainder) <= 0 && remainder == 0){
            printf("%d", remainder);
            return return_error("Argument error\nArgument [final_cluster_count] isn\'t in a proper format. \
It must be an integer higher than 0", EXIT_FAILURE);
        }
    }

    int len = load_clusters(argv[1], &clusters);

    //Checking if there was any error in load_clusters()
    if(len == 0){
        return EXIT_FAILURE;
    }

    //nearest neighbors algorithm implemented in the cycle below
    while(len > count){
        int c1, c2;
        find_neighbours(clusters, len, &c1, &c2);
        merge_clusters(&clusters[c1], &clusters[c2]);
        len = remove_cluster(clusters, len, c2);
    }
    print_clusters(clusters, len);

    free_clusters(&clusters, len);

    return 0;
}