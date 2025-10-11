import sys

MATCH =  +4
MISMATCH = -4
GAP = -8
NEG_INF = -10**12  # très grand négatif pour "interdire"

# lit les séquences depuis un fichier
def read_all_fastq(path):
    seqs = []
    with open(path, "r") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()   # ligne séquence
            f.readline()                 
            f.readline()                
            seqs.append(seq)
          
    return seqs

"""
    Chevauchement suffixe(x) vs préfixe(y) = semi-global asymétrique.
    Init : V(i,0)=0 (on peut sauter le préfixe de x gratuitement)
           V(0,j)=-inf pour j>0 (on DOIT commencer au début de y)
    Récurrence : max(diag, haut, gauche) avec scores +4/-4/-8
    Fin : score = max_j V(n, j) ; traceback à partir de (n, j*)
    """
def align_overlap(x, y):
    n = len(x)
    m = len(y)

    # Matrice des scores
    V = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    P = [[None] * (m + 1) for _ in range(n + 1)]  # pour le traceback

    # Initialisation
    V[0][0] = 0
    for i in range(1, n + 1):
        V[i][0] = 0           # départ libre dans x (préfixe ignoré)
        P[i][0] = 'STOP'      # marquer la bordure gauche comme point d'arrêt
    for j in range(1, m + 1):
        V[0][j] = NEG_INF     # interdit de sauter le début de y

    # Remplissage de la matrice
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = MATCH if x[i-1] == y[j-1] else MISMATCH
            diag = V[i-1][j-1] + score
            up = V[i-1][j] + GAP
            left = V[i][j-1] + GAP

            # Choix du max
            if diag >= up and diag >= left:
                V[i][j], P[i][j] = diag, 'DIAG'
            elif up >= left:
                V[i][j], P[i][j] = up, 'UP'
            else:
                V[i][j], P[i][j] = left, 'LEFT'
    
    # Meilleur score dans la dernière ligne
    j_star = max(range(m + 1), key=lambda j: V[n][j])
    best_score = V[n][j_star]

    # Traceback
    i, j = n, j_star
    a1, a2 = [], []
    while i > 0 and j >= 0:
        if j == 0 or i == 0 or P[i][j] == 'STOP':  # atteint le début de y → on a un préfixe(y)
            break
        move = P[i][j]
        if move == 'DIAG':
            a1.append(x[i-1]); a2.append(y[j-1])
            i -= 1; j -= 1
        elif move == 'UP':
            a1.append(x[i-1]); a2.append('-')
            i -= 1
        elif move == 'LEFT':
            a1.append('-');    a2.append(y[j-1])
            j -= 1
        else:  # sécurité
            break

    # inverser les alignements
    #align1 = ''.join(reversed(a1))
    #align2 = ''.join(reversed(a2))

    # Longueur du chevauchement = colonnes où les deux lettres sont réelles
    #overlap_len = sum(1 for c1, c2 in zip(align1, align2) if c1 != '-' and c2 != '-')

    return best_score

# Construire la matrice de chevauchement pour une liste de reads
def build_overlap_matrix(reads):
    n = len(reads)
    matrix = [[0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        for j in range(n):
            if i == j:
                matrix[i][j] = 0  # diagonale nulle
            else:
                score = align_overlap(reads[i], reads[j])
                matrix[i][j] = score
    return matrix
 
# Construire une liste d'arêtes (i, j, score) pour les scores >= threshold
def build_graph(matrix, threshold=80):
    n = len(matrix)
    edges = []
    for i in range(n):
        for j in range(n):
            if i != j and matrix[i][j] >= threshold:
                edges.append((i + 1, j + 1, matrix[i][j]))  # R1 -> R2 ...
    return edges

# exporter la liste d'arêtes au format DOT
def export_to_dot(edges, path="overlap_graph.dot"):
    with open(path, "w") as f:
        f.write("digraph G {\n")
        for i, j, score in edges:
            f.write(f'  R{i} -> R{j} [label="{score}"];\n')
        f.write("}\n")
    print(f"Graph exported to {path}")

def transitive_reduction(edges):
    #Construire une structure d'adjacence
    adj = {}
    for i, j, score in edges:
        adj.setdefault(i, []).append((j, score))

    # Fonction DFS pour vérifier si un chemin existe de src à dst
    def reachable(src, dst, visited=None):
        if visited is None:
            visited = set()
        if src == dst:
            return True
        visited.add(src)
        for (nxt, _) in adj.get(src, []):
            if nxt not in visited and reachable(nxt, dst, visited):
                return True
        return False

    reduced = []
    for (i, j, score) in edges:
        # Vérifie s’il existe un chemin indirect i -> k -> j
        keep = True
        for (k, _) in adj.get(i, []):
            if k != j and reachable(k, j):
                keep = False
                break
        if keep:
            reduced.append((i, j, score))

   

    return reduced

# Supprimer les cycles à deux nœuds en gardant toujours la meilleure direction
def collapse_two_cycles(edges):
    best = {}  # clé (min(i,j), max(i,j)) -> (src, dst, score)
    for i, j, s in edges:
        key = (min(i, j), max(i, j))
        if key not in best:
            best[key] = (i, j, s)
        else:
            i2, j2, s2 = best[key]
            # garder la direction avec le plus grand score
            if s > s2:
                best[key] = (i, j, s)
    # reconstruire la liste sans doublons
    final_edges = []
    seen = set()
    for i, j, s in best.values():
        if (i, j) not in seen:
            final_edges.append((i, j, s))
            seen.add((i, j))
    return final_edges


def assemble_sequence(reads, order, matrix):
    # Assemble la séquence finale à partir de l'ordre des reads et de la matrice des chevauchements
    assembled = reads[order[0] - 1]
    total_length = len(assembled)

    for i in range(len(order) - 1):
        a = order[i] - 1
        b = order[i + 1] - 1

        # recalcul du chevauchement réel entre read[a] et read[b]
        overlap_len = estimate_overlap(reads[a], reads[b])
        assembled += reads[b][overlap_len:]
        total_length = len(assembled)

    return assembled, total_length


def estimate_overlap(x, y):
    # Estime la longueur du chevauchement suffixe(x)/préfixe(y)
    max_len = min(len(x), len(y))
    for k in range(max_len, 0, -1):
        if x[-k:] == y[:k]:
            return k
    return 0


def main():
    path = "reads.fq"  # fichier FASTQ
    reads = read_all_fastq(path)
    n = len(reads)

    print("Nombre de reads lus :", n)
    print("Calcul de la matrice des chevauchements...\n")

    matrix = build_overlap_matrix(reads)

    print("MATRICE DES SCORES (20x20):")
    for i in range(n):
        for j in range(n):
            print(f"{matrix[i][j]:6}", end=" ")
        print()

    # Graphe filtré score >= 80
    edges = build_graph(matrix, threshold=80)
    print(f"\nNombre d'arêtes avec score >= 80 : {len(edges)}")
    export_to_dot(edges, path="overlap_graph.dot")

    # Réduction transitive
    edges = collapse_two_cycles(edges)
    reduced_edges = transitive_reduction(edges)
    

    print(f"Après réduction transitive et suppression des 2-cycles, nombre d'arêtes : {len(reduced_edges)}")
    export_to_dot(reduced_edges, path="overlap_graph_reduced.dot")

    print("Arêtes conservées (ordre des reads et chevauchements) :\n")
    for i, j, score in reduced_edges:
        print(f"R{i} -> R{j} [score={score}]")

    order = [2, 14, 11, 8, 3, 5, 1, 10, 18, 20, 15, 19, 13, 12, 16, 4, 9, 7, 6, 17]
    assembled_seq, total_length = assemble_sequence(reads, order, matrix)

    print("\nSéquence génomique reconstituée :")
    print(assembled_seq)
    print(f"\nLongueur totale du fragment : {total_length} bases")

if __name__ == "__main__":
    main()

