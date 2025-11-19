import math

"""
    1.1 - Lit un fichier FASTA et retourne une liste de tuples (id, sequence).
    
    db_path : chemin vers la banque de données FASTA (tRNAs.fasta)
   
    """
def load_fasta(db_path):
    sequences = []      # Liste finale des tuples (id, sequence)
    current_id = None   # Identifiant de la séquence courrante
    current_seq = []    # Liste des lignes de la séquence courrante

    # Ouvre le fichier FASTA en lecture
    with open(db_path, 'r') as file:
        for line in file:

            # Supprime les espaces blancs en début et fin de ligne
            line = line.strip()  
            if not line:
                continue  # Ignore les lignes vides

            # Détecte une nouvelle séquence
            if line.startswith('>'):
                # Si une seq precédente existe, la sauvegarde
                if current_id is not None:
                    sequences.append((current_id, ''.join(current_seq)))

                # On commence une nouvelle séquence
                current_id = line[1:]  # Enlève le '>' 
                current_seq = []       # Réinitialise la séquence courante

            else:
                # Ajoute la ligne à la séquence courante
                current_seq.append(line)

        # Sauvegarde la dernière séquence après la fin du fichier
        if current_id is not None:
            sequences.append((current_id, ''.join(current_seq)))
    return sequences

"""
    1.2 - Extrait tous les k-mers de la séquence qui respectent la graine.
    La graine peut contenir :
        '1' : position obligatoire (doit exister dans le k-mer)
        '0' : don't care (peu importe le caractère)

    sequence : séquence nucléotidique (string)
    seed     : string représentant la graine (ex: '11101101001')

    """
def extract_kmers(sequence, seed):

    k = len(seed)               # Longueur du k-mer 
    kmers = []                  # Liste des k-mers extraits

    # Parcourt la séquence pour extraire les k-mers
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]  # Extrait le k-mer courant
        
        # Vérifie si le k-mer respecte la graine
        valid = True              
        for j in range(k):
            if seed[j] == '1':
                # vérifie si le caractère est valide
                if kmer[j] not in 'ACGT':
                    valid = False
                    break
            # '0' : don't care, on ne fait rien
                
        if valid:
            kmers.append((i, kmer))    # Ajoute le k-mer valide à la liste

    return kmers

"""
    1.3 - Trouve les HSPs initiaux en cherchant chaque k-mer (de l'input)
    dans chaque séquence de la banque.

    kmers : liste de tuples (input_pos, kmer)
    database_sequences : liste de tuples (seq_id, seq)

    """
def find_hsps(kmers, database_sequences):
    hsps = []  # Liste des HSPs trouvés

    # Parcourt chaque k-mer
    for seq_id, seq in database_sequences:

        # Parcourt chaque k-mer de l'input
        for input_pos, kmer in kmers:
            k = len(kmer)

            # Cherche toutes les occurrences du k-mer dans la séquence de la banque
            for db_pos in range(len(seq) - k + 1):
                if seq[db_pos:db_pos + k] == kmer:
                    # Enregistre le HSP trouvé
                    hsps.append({
                        'seq_id': seq_id,
                        'input_pos': input_pos,
                        'db_pos': db_pos,
                        'kmer': kmer
                    })
    return hsps

"""
    1.4 - Étend un HSP dans les deux directions de manière gloutonne.
    
    input_seq : séquence input (string)
    db_seq    : séquence de la banque (string)
    input_start : début du HSP dans l'input
    db_start    : début du HSP dans la séquence de la banque
    k : taille du k-mer (graine)
    E : seuil (par défaut 4)
    
    """
def extend_hsp(input_seq, db_seq, input_start, db_start, k, E=4):

    # 1 - Initialisation du HSP
    input_left = input_start
    input_right = input_start + k - 1
    db_left = db_start
    db_right = db_start + k - 1

    # Score initial (k matches
    score = k * 5
    max_score = score  # Verifie si chute > E

    # 2 - Fonction pour comparer deux nucléotides
    def score_char(a,b):
        return 5 if a == b else -4
    
    # 3 - Extension gloutonne 
    while True:
        # Veifie si on peut étendre à gauche
        can_extend_left = (input_left > 0) and (db_left > 0)

        # Veifie si on peut étendre à droite
        can_extend_right = (input_right < len(input_seq) - 1) and (db_right < len(db_seq) - 1)

        # Si on ne peut pas étendre, on arrête
        if not can_extend_left and not can_extend_right:
            break

        # Calcul potentiel score gauche
        left_score = None
        if can_extend_left:
            a = input_seq[input_left - 1]
            b = db_seq[db_left - 1]
            left_score = score_char(a,b)

        # Calcul potentiel score droite
        right_score = None
        if can_extend_right:
            a = input_seq[input_right + 1]
            b = db_seq[db_right + 1]
            right_score = score_char(a,b)

        # Décide l'extension
        # Cas: les deux possibles
        if left_score is not None and right_score is not None:
            if left_score > right_score:
                direction = "L"
            elif right_score > left_score:
                direction = "R"
            else:
                direction = "L"  

        # Cas : seulement gauche possible
        elif left_score is not None:
            direction = "L"
        # Cas : seulement droite possible
        else:
            direction = "R"

        # Effectue l'extension choisie
        if direction == "L":
            new_score = score + left_score
            new_input_left = input_left - 1
            new_db_left = db_left - 1

            # Vérifie du seuil -E
            if max_score - new_score > E:
                # On arrête complètement l'extension
                break

            # Validation de l'extension
            input_left = new_input_left
            db_left = new_db_left
            score = new_score

        else:  # direction == "R"
            new_score = score + right_score
            new_input_right = input_right + 1
            new_db_right = db_right + 1

            # Vérifie du seuil -E
            if max_score - new_score > E:
                # On arrête complètement l'extension
                break

            # Validation de l'extension
            input_right = new_input_right
            db_right = new_db_right
            score = new_score

        # Met à jour le score max
        if score > max_score:
            max_score = score

    # 4 - Retourne les résultats de l'extension
    align_input = input_seq[input_left:input_right + 1]
    align_db = db_seq[db_left:db_right + 1]

    # 5 - Retourne les résultats
    return {
        'input_start': input_left,
        'input_end': input_right,
        'db_start': db_left,
        'db_end': db_right,
        'score': score,
        'align_input': align_input,
        'align_db': align_db
    }

"""
    1.5 - Fusionne les HSPs chevauchants pour chaque séquence de la banque.

    hsps : liste de dictionnaires HSPs après extension :
        {
            'seq_id': ...,
            'input_start': ...,
            'input_end': ...,
            'db_start': ...,
            'db_end': ...,
            'score': ...,
            'align_input': ...,
            'align_db': ...
        }
    """
def merge_hsps(hsps):
    # 1 - Regroupe les HSPs par seq_id
    from collections import defaultdict
    grouped = defaultdict(list)

    for h in hsps:
        grouped[h['seq_id']].append(h)
    
    merged_hsps = []

    # 2 - Pour chaque groupe
    for seq_id, hsp_list in grouped.items():

        # Trie par position db_start
        hsp_list.sort(key=lambda x: x['db_start'])

        merged = []
        current = hsp_list[0]  # Hsp courant à fusionner

        # 3 - Fusionne les HSPs
        for h in hsp_list[1:]:
            # Vérifie le chevauchement
            if h['db_start'] <= current['db_end']:
                # Fusionne
                current['db_end'] = max(current['db_end'], h['db_end'])
                current['input_start'] = min(current['input_start'], h['input_start'])
                current['input_end'] = max(current['input_end'], h['input_end'])
                current['score'] = max(current['score'], h['score'])  

            else:
                # Pas de chevauchement, sauvegarde le courant
                merged.append(current)
                current = h  # Nouveau courant

        # Sauvegarde le dernier HSP
        merged.append(current)

        # Ajoute au resultat final
        merged_hsps.extend(merged)
    return merged_hsps
        
"""
    1.6 - Calcule bitscore + e-value pour chaque HSP, puis filtre les HSPs significatifs.
    
    hsps : liste d'HSPs après extension et fusion
    total_db_size : m = taille totale des séquences de la banque
    input_length  : n = taille de la séquence input
    ss_threshold  : seuil -ss (par défaut 1e-3)

    """
def score_and_filter_hsps(hsps, total_db_size, input_length, ss_threshold=1e-3):
    λ = 0.192
    K = 0.176

    # Ajouter bitscore et e-value à chaque HSP
    for h in hsps:
        S = h['score']

        # 1 - Calcul du bitscore
        B = (λ * S - math.log(K)) / math.log(2)
        h['bitscore'] = round(B)

        # 2 - Calcul de l'e-value
        # e = m *n * 2^(-B)
        e_value = total_db_size * input_length * (2 ** (-h['bitscore']))
        h['evalue'] = e_value

    # 3 - Filtrage des HSPs significatifs
    significant = [h for h in hsps if h['evalue'] <= ss_threshold]

    # S"il y en aucun 
    if not significant:
        return []

    # 4 - Garder un seul HSP par sequence
    from collections import defaultdict
    best_per_seq = defaultdict(list)

    for h in significant:
        best_per_seq[h['seq_id']].append(h)

    final_hsps = []

    for seq_id, hsp_list in best_per_seq.items():
        # Trie par bitscore décroissant
        hsp_list.sort(key=lambda x: x['bitscore'], reverse=True)
        # Garde le meilleur
        final_hsps.append(hsp_list[0])

    return final_hsps

"""
1.7 – Produit l’output final : sélectionne et affiche, pour chaque séquence
de la banque, le meilleur HSP significatif (trié par pertinence).
"""
def plast_search(input_seq, db_path, seed="11111111111", E=4, ss_threshold=1e-3):
    # 1 - Charge la banque de données FASTA
    database = load_fasta(db_path)

    # Taille totale des séquences de la banque
    total_db_size = sum(len(seq) for _, seq in database)
    input_length = len(input_seq)

    # 2 - Extrait les k-mers de la séquence input
    kmers = extract_kmers(input_seq, seed)

    # 3 - Trouve les HSPs initiaux
    initial_hsps = find_hsps(kmers, database)

    # 4 - Étend chaque HSP
    extended_hsps = []
    for h in initial_hsps:
        seq_id = h['seq_id']
        kmer = h['kmer']
        input_start = h['input_pos']
        db_start = h['db_pos']

        # On recupere la sequence correspondante dans la banque
        db_seq = next(seq for (id_, seq) in database if id_ == seq_id)

        extended = extend_hsp(
            input_seq,
            db_seq,
            input_start,
            db_start,
            k=len(kmer),
            E=E
        )

        extended['seq_id'] = seq_id  # Ajoute l'id de la séquence
        extended_hsps.append(extended)


    # 5 - Fusionne les HSPs chevauchants
    merged = merge_hsps(extended_hsps)

    # 6 - Calcule bitscore + e-value et filtre les HSPs significatifs
    final_hsps = score_and_filter_hsps(
        merged,
        total_db_size,
        input_length,
        ss_threshold
    )

    # 7 - Trie les HSPs finaux par bitscore décroissant
    final_hsps.sort(key=lambda x: x['bitscore'], reverse=True)

    # 8 - Affiche les résultats
    print("\n -- Resultats PLAST --")
    if not final_hsps:
        print("Aucun HSP significatif trouvé (e-value > seuil).")
        return
    
    for h in final_hsps:
        print("\nSequence :", h['seq_id'])
        print("Bitscore :", h['bitscore'])
        print("E-value  :", h['evalue'])
        print("Input    : [{}:{}]".format(h['input_start'], h['input_end']))
        print("DB       : [{}:{}]".format(h['db_start'], h['db_end']))
        print("Align input :", h['align_input'])
        print("Align db    :", h['align_db'])

# -- Exécution principale --
if __name__ == "__main__":

    #  1) Charger la séquence à rechercher 
    input_path = "unknown.fasta"
    input_data = load_fasta(input_path)
    input_seq = input_data[0][1]  # On prend la première séquence du fichier

    # 2) Charger la banque de données 
    db_path = "tRNAs.fasta"

    # 3) Paramètres simples
    seed = "11111111111"   # Graine par défaut
    E = 4                  # Seuil d’extension
    ss_threshold = 1e-3    # Seuil de significativité

    # 4) Exécuter le PLAST 
    plast_search(
        input_seq,
        db_path=db_path,
        seed=seed,
        E=E,
        ss_threshold=ss_threshold
    )
