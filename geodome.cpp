//
// grantham@plunk.org, April 2011
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if _WIN32
#define M_PI 3.141592
#endif

#undef OCTAHEDRON

int triangles[][3] =
{
#ifdef OCTAHEDRON
    {0, 1, 2},
    {0, 2, 3},
    {0, 3, 4},
    {0, 4, 1},
    {5, 2, 1},
    {5, 3, 2},
    {5, 4, 3},
    {5, 1, 4},
#else
    {0, 1, 2},
    {0, 2, 3},
    {0, 3, 4},
    {0, 4, 5},
    {0, 5, 1},
    {1, 6, 2},
    {6, 7, 2},
    {2, 7, 3},
    {7, 8, 3},
    {3, 8, 4},
    {8, 9, 4},
    {4, 9, 5},
    {9, 10, 5},
    {5, 10, 1},
    {10, 6, 1},
    {6, 11, 7},
    {7, 11, 8},
    {8, 11, 9},
    {9, 11, 10},
    {10, 11, 6},
#endif
};
const int src_triangle_count = sizeof(triangles) / sizeof(triangles[0]);
int triangle_count;

struct vertex {
    float m_v[3];
    float operator[](int i) {return m_v[i];};
    vertex operator-(const vertex &v) const
    {
        return vertex(m_v[0] - v.m_v[0], m_v[1] - v.m_v[1], m_v[2] - v.m_v[2]);
    }
    vertex operator+(const vertex &v) const
    {
        return vertex(m_v[0] + v.m_v[0], m_v[1] + v.m_v[1], m_v[2] + v.m_v[2]);
    }
    vertex operator*(float f) const
    {
        return vertex(m_v[0] * f, m_v[1] * f, m_v[2] * f);
    }
    vertex operator/(float f) const
    {
        return vertex(m_v[0] / f, m_v[1] / f, m_v[2] / f);
    }
    vertex() {};
    vertex(float x, float y, float z)
    {
        m_v[0] = x; m_v[1] = y; m_v[2] = z;
    }
    vertex &normalize()
    {
        float d;
        d = sqrt(m_v[0] * m_v[0] + m_v[1] * m_v[1] + m_v[2] * m_v[2]);
        m_v[0] /= d;
        m_v[1] /= d;
        m_v[2] /= d;
        return *this;
    }
    float length() const {
        return sqrt(m_v[0] * m_v[0] + m_v[1] * m_v[1] + m_v[2] * m_v[2]);
    }
    bool closeto(const vertex &v, float tolerance = .0001) const
    {
        vertex d = *this - v;
        float dist = d.length();
        return dist < tolerance;
    }
    bool closeto(vertex &v, float tolerance = .0001) const
    {
        vertex d = *this - v;
        float dist = d.length();
        return dist < tolerance;
    }
} vertices[12];

void emit_vertex(vertex &v)
{
    printf("%f %f %f", v[0], v[1], v[2]);
}

struct triangle {
    vertex m_v[3];
    vertex &operator[](int i) { return m_v[i]; }
    const vertex &operator[](int i) const { return m_v[i]; }
};

vertex make_flat(vertex n0, vertex n1, vertex n2)
{
    int i;
    float d = 0;
    vertex n;

    n = (n0 + n1 + n2) / 3.0;
    return n / n.length();
}

vertex get_barycenter(const vertex &v0, const vertex &v1, const vertex &v2)
{
    return (v0 + v1 + v2) / 3.0;
}

void emit_readable(triangle *tris, int count)
{
    int i, j, k;
    vertex n;
    float a;
    vertex vertices[count * 3];
    int vertex_count = 0;
    int edges[count * 3][2];
    int edge_count = 0;
    int triangles[count][3];

    for(i = 0; i < count; i++) {
        triangle &t = tris[i];
	int elements[3];

	for(j = 0; j < 3; j++) {
	    for(k = 0; k < vertex_count; k++)
		if(tris[i][j].closeto(vertices[k], .001))
		    break;
	    if(k == vertex_count) {
		vertices[k] = tris[i][j];
		vertex_count++;
	    }
	    elements[j] = k;
	    triangles[i][j] = elements[j];
	}
	for(j = 0; j < 3; j++) {
	    int i0, i1;
	    i0 = elements[j];
	    i1 = elements[(j + 1) % 3];
	    for(k = 0; k < edge_count; k++)
		if( (edges[k][0] == i0 && edges[k][1] == i1) ||
		    (edges[k][0] == i1 && edges[k][1] == i0) )
		    break;
	    if(k == edge_count) {
		edges[k][0] = i0;
		edges[k][1] = i1;
		edge_count++;
	    }
	}
    }

    printf("# A vertex is the same as a \"corner\"\n");
    printf("vertices %d\n", vertex_count);
    for(i = 0; i < vertex_count; i++) {
	printf("vertex ");
	emit_vertex(vertices[i]);
	puts("");
    }
    printf("\n# Edges, the 0-based indices of the vertices at each end\n");
    printf("edges %d\n", edge_count);
    for(i = 0; i < edge_count; i++) {
	printf("edge %d %d\n", edges[i][0], edges[i][1]);
    }
    printf("\n# Triangles, the 0-based indices of the vertices\n");
    printf("triangles %d\n", count);
    for(i = 0; i < count; i++) {
	printf("triangle %d %d %d\n", triangles[i][0], triangles[i][1],
	    triangles[i][2]);
    }
}


void emit_triangles(triangle *tris, int count)
{
    int i;
    vertex n;
    float a;

    for(i = 0; i < count; i++) {
        triangle &t = tris[i];

        n = make_flat(t[0], t[1], t[2]);

        a = get_barycenter(t[0], t[1], t[2])[1] * .5 + .5;

        printf("\"*\" none 0 0 0 1 100 ");

        emit_vertex(t[0]); printf(" ");
        emit_vertex(n); printf(" ");
        // printf("1 1 1 1 0 0 ");
        printf("0 0 0 1 0 0 ");

        emit_vertex(t[1]); printf(" ");
        emit_vertex(n); printf(" ");
        // printf("1 1 1 1 0 0 ");
        printf("0 0 0 1 0 0 ");

        emit_vertex(t[2]); printf(" ");
        emit_vertex(n); printf(" ");
        // printf("1 1 1 1 0 0 ");
        printf("0 0 0 1 0 0 ");
        puts("");
    }
}

bool do_normalize = true;

void dice_triangle(triangle &in, int level, triangle* out)
{
    int i, j;
    vertex v01inc;
    vertex v12inc;
    triangle* t = out;

    v01inc = (in[1] - in[0]) / (float)level;
    v12inc = (in[2] - in[1]) / (float)level;

    for(i = 0; i < level; i++) {
        vertex top = in[0] + v01inc * i;
        vertex bot = in[0] + v01inc * (i + 1);
        // construct swath parameters
        for(j = 0; j <= i; j++) {
            // iterate down swath
            t->m_v[0] = top + v12inc * j;
            if(do_normalize) t->m_v[0].normalize();
            t->m_v[1] = bot + v12inc * j;
            if(do_normalize) t->m_v[1].normalize();
            t->m_v[2] = bot + v12inc * (j + 1);
            if(do_normalize) t->m_v[2].normalize();

            t++;

            // skip the last one
            if(j == i)
                break;

            t->m_v[0] = top + v12inc * j;
            if(do_normalize) t->m_v[0].normalize();
            t->m_v[1] = bot + v12inc * (j + 1);
            if(do_normalize) t->m_v[1].normalize();
            t->m_v[2] = top + v12inc * (j + 1);
            if(do_normalize) t->m_v[2].normalize();
            t++;
        }
    }
}

void tessellate_triangles(triangle *triangles_in, triangle **triangles_out,
    int triangle_count, int *new_count,
    int level)
{
    int i;

    *new_count = triangle_count * level * level;
    *triangles_out = (triangle*)malloc(sizeof(triangle) * *new_count);

    for(i = 0; i < triangle_count; i++)
        dice_triangle(triangles_in[i], level,
            (*triangles_out) + i * level * level);
}

struct edge {
    vertex v[2];
    vertex &operator[](int i) { return v[i]; }
    const vertex &operator[](int i) const { return v[i]; }
    bool operator ==(const edge &e) const {
        return (v[0].closeto(e[0]) && v[1].closeto(e[1])) ||
            (v[0].closeto(e[1]) && v[1].closeto(e[0]));
    }
    bool operator ==(edge &e) {
        return (v[0].closeto(e[0]) && v[1].closeto(e[1])) ||
            (v[0].closeto(e[1]) && v[1].closeto(e[0]));
    }
    void set(vertex &v1, vertex &v2) { v[0] = v1; v[1] = v2; }
    edge(vertex &v1, vertex &v2) { set(v1, v2); }
};

void emit_edges(edge *edges, int count, char *tag)
{
    int i;
    vertex n;
    float a;

    printf("int %sCount = %d;\n", tag, count);
    printf("float %sFloats[][2][3] = {\n", tag);

    for(i = 0; i < count; i++) {
        edge &e = edges[i];
        vertex &v0 = e[0];
        vertex &v1 = e[1];
        printf("    { {%f, %f, %f}, {%f, %f, %f} },\n",
            v0[0], v0[1], v0[2],
            v1[0], v1[1], v1[2]);
    }
    printf("};\n", tag);
}

int sort_two_triangles(const void *p1, const void *p2)
{
    const triangle &t1 = *(const triangle *)p1;
    const triangle &t2 = *(const triangle *)p2;

    vertex c1 = get_barycenter(t1[0], t1[1], t1[2]);
    vertex c2 = get_barycenter(t2[0], t2[1], t2[2]);

    if(c1[1] < c2[1])
        return 1;
    if(c1[1] == c2[1])
        return 0;
    return -1;
}

void sort_barycenters(triangle *tris, int count)
{
    qsort(tris, count, sizeof(triangle), sort_two_triangles);
}

int sort_two_edges(const void *p1, const void *p2)
{
    const edge &e1 = *(const edge *)p1;
    const edge &e2 = *(const edge *)p2;

    vertex c1 = (e1[0] + e1[1]) / 2.0;
    vertex c2 = (e2[0] + e2[1]) / 2.0;

    if(c1[1] < c2[1])
        return 1;
    if(c1[1] == c2[1])
        return 0;
    return -1;
}

void sort_edge_centers(edge *edges, int count)
{
    qsort(edges, count, sizeof(edge), sort_two_edges);
}

void make_octahedron_verts()
{
    int i = 0;

    vertices[i++] = vertex(0, 1, 0);
    vertices[i++] = vertex(-1, 0, 0);
    vertices[i++] = vertex(0, 0, 1);
    vertices[i++] = vertex(1, 0, 0);
    vertices[i++] = vertex(0, 0, -1);
    vertices[i++] = vertex(0, -1, 0);
}

void make_icosahedron_verts()
{
    int i;

    vertices[0] = vertex(0, 1, 0);

    for(i = 0; i < 5; i++) {
        float alpha, beta, x, y, z;
        alpha = i * (2 * M_PI) / 5;
        beta = 1.1072;
        y = cos(beta);
        x = cos(alpha) * sin(beta);
        z = sin(alpha) * sin(beta);
        vertices[i + 1] = vertex(x, y, z);
    }
    for(i = 0; i < 5; i++) {
        float alpha, beta, x, y, z;
        alpha = (i + 0.5) * (2 * M_PI) / 5;
        beta = M_PI - 1.1072;
        y = cos(beta);
        x = cos(alpha) * sin(beta);
        z = sin(alpha) * sin(beta);
        vertices[i + 6] = vertex(x, y, z);
    }
    vertices[11] = vertex(0, -1, 0);
}

void make_edges(triangle *tris, int triangle_count,
    edge **edgesReturn, int *edge_countReturn)
{
    int i, j, k;
    edge *edges;
    int edge_count = 0;

    edges = (edge*)malloc(sizeof(edge) * triangle_count * 3);

    for(i = 0; i < triangle_count; i++) {
        triangle &t = tris[i];
        for(j = 0; j < 3; j++) {
            int v1, v2;
            v1 = j;
            v2 = (j + 1) % 3;
            edge e(t[v1], t[v2]);
            for(k = 0; k < edge_count; k++)
                if(e == edges[k])
                    break;
            if(k == edge_count)
            {
                edges[k] = e;
                edge_count++;
            }
        }
    }
    *edgesReturn = edges;
    *edge_countReturn = edge_count;
}

void usage(char *progname)
{
    fprintf(stderr, "Usage: %s [options] tesselation-factor\n", progname);
    fprintf(stderr, "Options\n");
    fprintf(stderr, "\t--trisrc        Output \"trisrc\" format data instead\n");
    fprintf(stderr, "\t                of text file format\n");
    fprintf(stderr, "\t--no-normalize  Don't normalize resulting triangles\n");
    fprintf(stderr, "\t--original N    Output only triangles resulting from original\n");
    fprintf(stderr, "\t                triangle N\n");
    fprintf(stderr, "\t--top N         Output only triangles resulting from first N triangles\n");
}

int main(int argc, char **argv)
{
    int j;
    int i;
    int level = 3;
    bool only_one_source_tri = false;
    int which_triangle;
    int tri_limit = -1;
    bool emit_text_file = true;
    char *progname = argv[0];

    argc--;
    argv++;

    while(argc > 0 && argv[0][0] == '-') {
	if(strcmp(argv[0], "--no-normalize") == 0) {
	    do_normalize = false;
	    argv++;
	    argc--;
	} else if(strcmp(argv[0], "--trisrc") == 0) {
	    emit_text_file = false;
	    argv++;
	    argc--;
	} else if(strcmp(argv[0], "--original") == 0) {
	    if(argc < 2) {
		fprintf(stderr, "expected triangle number for --original\n");
                usage(progname);
		exit(EXIT_FAILURE);
	    }
	    only_one_source_tri = true;
	    which_triangle = atoi(argv[1]);
	    argv += 2;
	    argc -= 2;
	} else if(strcmp(argv[0], "--top") == 0) {
	    if(argc < 2) {
		fprintf(stderr, "expected triangle number for --top\n");
                usage(progname);
		exit(EXIT_FAILURE);
	    }
	    tri_limit = atoi(argv[1]);
	    argv += 2;
	    argc -= 2;
	} else {
	    fprintf(stderr, "don't know what to do with \"%s\"\n", argv[0]);
            usage(progname);
	    exit(EXIT_FAILURE);
	}
    }

    if(argc < 1) {
	fprintf(stderr, "expected subdivision factor\n");
        usage(progname);
	exit(EXIT_FAILURE);
    }

    level = atoi(argv[0]);

#ifdef OCTAHEDRON
    make_octahedron_verts();
#else
    make_icosahedron_verts();
#endif

    triangle *ico = (triangle *)malloc(sizeof(triangle) * src_triangle_count);
    triangle_count = 0;
    for(i = 0; i < src_triangle_count; i++) {
	if(only_one_source_tri && i != which_triangle)
	    continue;
        for(j = 0; j < 3; j++)
            ico[triangle_count][j] = vertices[triangles[i][j]];
	triangle_count++;
    }

    triangle *tessellated;
    int tessellated_count;
    tessellate_triangles(ico, &tessellated, triangle_count, &tessellated_count, level);

    sort_barycenters(tessellated, tessellated_count);

    if(emit_text_file)
	emit_readable(tessellated, (tri_limit > -1) ? tri_limit : tessellated_count);
    else
	emit_triangles(tessellated, (tri_limit > -1) ? tri_limit : tessellated_count);

    free(tessellated);
    free(ico);
}
