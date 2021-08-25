
(* Sequences genetiques *)

type nucleotide = A | C | G | T;;
type brin = nucleotide list;;



(** Question 1 **)

(* Fonction retournant le contenu GC d'un brin d'ADN passé en argument *)
let contenu_gc (x : brin) : float = List.fold_left (
    fun compteur e -> if (fun e -> e = C || e = G) e 
      then compteur +. (1./. float_of_int (List.length x)) 
      else compteur
  ) 0. x;;

(* TEST *)
let _ = assert (contenu_gc [A; T; G; T; T; G; A; C] = 0.375);;
let _ = assert (contenu_gc [C; T; T; A] = 0.25);;
let _ = assert (contenu_gc [A; A; A; T; A] = 0.);;



(** Question 2 **)
(* Fonction calculant et retournant le brin complémentaire du brin donné en paramètre *)
let rec brin_complementaire (x : brin) : brin = match x with
  | [] -> []
  | x::r -> if x = A 
    then T::brin_complementaire r 
    else 
    if x = T 
    then A::brin_complementaire r 
    else
    if x = C 
    then G::brin_complementaire r 
    else 
    if x = G then C::brin_complementaire r else [C];;

(* TEST *)
let _ = assert (brin_complementaire [T] = [A]);;
let _ = assert (brin_complementaire [C; T; T; C] = [G; A; A; G]);;
let _ = assert (brin_complementaire [C; T; A; A; T; G; T] = [G; A; T; T; A; C; A]);;



(** Question 3 **)

(* Fonction retournant la distance d'édition entre les 2 brins donnés en paramètres *)
let distance (x : brin) (y : brin) : int = 
  List.fold_left2 (fun acc x y -> if x = y then acc else 1 + acc) 0 x y;;

(* TEST *)
let _ = assert (distance [T] [T] = 0);;
let _ = assert (distance [T] [C] = 1);;
let _ = assert (distance [G; A; G] [A; G; G] = 2);;




(** Question 4 **)

(* Fonction calculant et retournant la similarité entre 2 brins donnés en paramètre *)
let similarite (x : brin) (y : brin) : float = 
  if (List.length x) <> (List.length y)
  then raise (Failure "brin invalide")
  else
    1.0 -. (float_of_int (distance x y)) /. (float_of_int (List.length x));;

(* TEST *)
let _ = assert (similarite [C; G; A; T] [T; A; G; T] = 0.25);;
let _ = assert (similarite [A; G; C; T] [T; A; A; G] = 0.);;
let _ = assert (similarite [A; G; C; T] [A; G; C; T] = 1.);;
let _ = assert (similarite [A; G; C; T] [A; G; C] = raise (Failure "brin invalide"));;


type acide = Ala | Arg | Asn | Asp | Cys
           | Glu | Gln | Gly | His | Ile
           | Leu | Lys | Phe | Pro | Ser
           | Thr | Trp | Tyr | Val | START | STOP

let codon_vers_acide (n1 : nucleotide) (n2 : nucleotide) (n3 : nucleotide) : acide = 
  match (n1, n2, n3) with
  | (A, A, A) -> Phe | (A, A, G) -> Phe | (A, A, T) -> Leu  | (A, A, C) -> Leu
  | (G, A, A) -> Leu | (G, A, G) -> Leu | (G, A, T) -> Leu  | (G, A, C) -> Leu
  | (T, A, A) -> Ile | (T, A, G) -> Ile | (T, A, T) -> Ile  | (T, A, C) -> START
  | (C, A, A) -> Val | (C, A, G) -> Val | (C, A, T) -> Val  | (C, A, C) -> Val
  | (A, G, A) -> Ser | (A, G, G) -> Ser | (A, G, T) -> Ser  | (A, G, C) -> Ser
  | (G, G, A) -> Pro | (G, G, G) -> Pro | (G, G, T) -> Pro  | (G, G, C) -> Pro
  | (T, G, A) -> Thr | (T, G, G) -> Thr | (T, G, T) -> Thr  | (T, G, C) -> Thr
  | (C, G, A) -> Ala | (C, G, G) -> Ala | (C, G, T) -> Ala  | (C, G, C) -> Ala
  | (A, T, A) -> Tyr | (A, T, G) -> Tyr | (A, T, T) -> STOP | (A, T, C) -> STOP
  | (G, T, A) -> His | (G, T, G) -> His | (G, T, T) -> Gln  | (G, T, C) -> Gln
  | (T, T, A) -> Asn | (T, T, G) -> Asn | (T, T, T) -> Lys  | (T, T, C) -> Lys
  | (C, T, A) -> Asp | (C, T, G) -> Asp | (C, T, T) -> Glu  | (C, T, C) -> Glu
  | (A, C, A) -> Cys | (A, C, G) -> Cys | (A, C, T) -> STOP | (A, C, C) -> Trp
  | (G, C, A) -> Arg | (G, C, G) -> Arg | (G, C, T) -> Arg  | (G, C, C) -> Arg
  | (T, C, A) -> Ser | (T, C, G) -> Ser | (T, C, T) -> Arg  | (T, C, C) -> Arg
  | (C, C, A) -> Gly | (C, C, G) -> Gly | (C, C, T) -> Gly  | (C, C, C) -> Gly;;

(** Question 5 **)

(* 
  Fonction auxiliaire vérifiant si la taille du brin den paramètre est bien un multiple de 3
  (pour une cohérence avec la taille des codons)
 *)
let verif_taille_brin (brin:brin): bool  = if List.length brin mod 3 = 0 then true else false;;

(* Fonction qui détermine et renvoie une liste d'acide aminés de la première chaîne du brin *)
let rec brin_vers_chaine (x : brin) : acide list = 
  if verif_taille_brin x = false 
  then raise (Failure "brin invalide")
  else match x with
    | e::y::z::l -> 
      if codon_vers_acide e y z = STOP 
      then [] 
      else 
      if codon_vers_acide e y z != START 
      then [codon_vers_acide e y z] @ brin_vers_chaine l 
      else brin_vers_chaine l
    | _ -> raise (Failure "brin invalide");;

(* TEST *)
let _ = assert (brin_vers_chaine [T; A; C; G; G; C; T; A; G; A; T; T; T; A; C; G; C; T; A; A; T; A; T; C] = [Pro; Ile]);;
let _ = assert (brin_vers_chaine [T; A; C; T; A; C] = raise (Failure "brin invalide"));;
let _ = assert (brin_vers_chaine [T; A; C; G; G; A; T; C] = raise (Failure "brin invalide"));;



(** Question 6 **)

(* Fonction retournant le décodage de toutes les chaînes du brin donné en paramètre *)
let rec brin_vers_chaines (x : brin) : acide list list  = 
  if verif_taille_brin x = false 
  then raise (Failure "brin invalide")
  else match x with
    | e::y::z::l -> if codon_vers_acide e y z = START then 
        [brin_vers_chaine l]@ brin_vers_chaines l 
      else brin_vers_chaines l
    | _ -> [];;

(* TEST *)
let _ = assert (brin_vers_chaines [T; A; C; G; G; C; T; A; G; A; T; T; T; A; C; G; C; T; A; A; T; A; T; C] = [[Pro; Ile]; [Arg; Leu]]);;
let _ = assert (brin_vers_chaines [T; A; C; T; A; C] = raise (Failure "brin invalide"));;
let _ = assert (brin_vers_chaines [T; A; C; G; G; A; T; C] = raise (Failure "brin invalide"));;

(* Arbres phylogenetiques *)

type arbre_phylo =
  | Lf of brin
  | Br of arbre_phylo * brin * int * arbre_phylo;;

(* Arbres de test *)
let arbre = Br(
    Br(Lf([G; C; A; T]), [A; C; A; T], 3, Lf([T; C; G; T])),
    [A; A; A; A],
    8,
    Br(Lf([T; A; G; A]), [A; A; G; A], 2, Lf([G; A; G; A]))
  );;

let arbre2 = Br(
    Br(Lf([G; A; A; T]), [G; C; T; T], 5, Lf([C; A; G; T])),
    [A; A; T; A],
    12,
    Br(Lf([T; A; G; A]), [A; A; G; A], 3, Lf([G; T; G; A]))
  );;

let arbre3 = Br(
    Br(Lf([G; C; A; T]), [A; C; A; T], 3, Lf([T; C; G; T])),
    [A; T; A; A],
    10,
    Br(Lf([T; A; G; A]), [A; A; G; A], 2, Lf([G; A; G; A]))
  );;
(** Question 1 **)

(* Fonction auxiliaire qui retourne une string formatée pour un brin en particulier *)
let rec affichage_brin (brin: brin) : string = match brin with
  | [] -> ""
  | x::s -> if x = A then "A"^(affichage_brin s) 
    else if x = T then "T"^(affichage_brin s) 
    else if x = C then "C"^(affichage_brin s) 
    else "G"^(affichage_brin s);;

(* Fonction retournant une string unique représentant de manière unique un arbre phylogénétique *)
let rec arbre_phylo_vers_string (a : arbre_phylo) : string = match a with 
  | Lf(brin) ->  (affichage_brin brin)
  | Br(l, brin, malus, r) ->  
    (arbre_phylo_vers_string l)^"\n"^(affichage_brin brin)^" ("^(string_of_int malus)^")"^"\n"^(arbre_phylo_vers_string r);;

(*** TEST ***)
let _ = assert (arbre_phylo_vers_string arbre = "GCAT\nACAT (3)\nTCGT\nAAAA (8)\nTAGA\nAAGA (2)\nGAGA");;
let _ = assert (arbre_phylo_vers_string arbre2 = "GAAT\nGCTT (5)\nCAGT\nAATA (12)\nTAGA\nAAGA (3)\nGTGA");;
print_string (arbre_phylo_vers_string arbre);;

(** Question 2 **)

(* Fonction de parcours d'un arbre (parcours racine - sous-arbre gauche - sous-arbre droit) *)
let rec arbre_parcours (arbre:arbre_phylo) : brin list = match arbre with
  | Lf(x) -> [x]
  | Br(l, x, malus,  r) ->  [x] @ arbre_parcours l @ arbre_parcours r ;;

(* Fonction calculant la similarité entre 2 arbres phylogénétiques *)
let arbre_similarite (a1:arbre_phylo) (a2:arbre_phylo) : float = 
  List.fold_left2 (fun acc x y -> acc +. similarite x y) 0. (arbre_parcours a1) (arbre_parcours a2);;

(* Fonction comparant les similarités de 2 arbres par rapport à un arbre de référence*)
let comparaison_similarite (ref:arbre_phylo) (a1:arbre_phylo) (a2:arbre_phylo) : arbre_phylo = 
  if abs_float (arbre_similarite ref a1) < abs_float (arbre_similarite ref a2) 
  then a1 else a2;;

(* Fonction retournant l'arbre le plus similaire à l'arbre a parmi une liste l d'arbres phylo *)
let rec similaire (a : arbre_phylo) (l : arbre_phylo list) : arbre_phylo = match l with
  | [] -> failwith "liste d'arbres invalide"
  | [arbre] -> arbre
  | x::suite -> comparaison_similarite a x (similaire a suite);;

(* TEST *)
let _ = assert (similaire arbre [arbre2; arbre] = arbre2);;


(** Questions 3 **)

(** Question 3.1 **)

(* Fonction d'extraction du brin de la racine de l'arbre *)
let get_root (a : arbre_phylo) : brin = match a with
  | Lf(x) -> x
  | Br(l, x, m, r) -> x;;

(* TEST *)
let _ = assert ((get_root arbre) = [A; A; A; A]);;
let _ = assert ((get_root arbre3) = [A; T; A; A]);;

(** Question 3.2 **)

(* Fonction d'extraction du malus de la racine de l'arbre *)
let get_malus (a : arbre_phylo) : int = match a with
  | Lf(x) -> 0
  | Br(l, x, m, r) -> m;;

(* TEST *)
let _ = assert ((get_malus arbre) = 8);;
let _ = assert ((get_malus arbre3) = 10);;

(** Question 3.3 **)

(* Fonction de construction d'un arbre phylo à partir de données en paramètre *)
let br (ag : arbre_phylo) (b : brin) (ad : arbre_phylo) : arbre_phylo = 
  Br(
    ag, 
    b, 
    distance (get_root ag) b + distance (get_root ad) b + (get_malus ag) + (get_malus ad), 
    ad);;

(* TEST *)
let _ = assert (br (Lf [A]) [T] (Lf [C]) = Br (Lf [A], [T], 2, Lf [C]));;
let _ = assert (br (br (Lf [T]) [G] (Lf [A])) [T] (br (Lf [T]) [G] (Lf [A])) = Br (Br (Lf [T], [G], 2, Lf [A]), [T], 6, Br (Lf [T], [G], 2, Lf [A])));;


(** Question 4 **)

(* Fonction retournant l'arbre avec le plus petit malus global parmi 2 arbres *)
let min_malus_comparaison (a1:arbre_phylo) (a2:arbre_phylo) : arbre_phylo =
  if get_malus a1 < get_malus a2 then a1 else a2;;

(* Fonction retournant l'arbre avec le plus petit malus global parmi une liste d'arbres *)
let rec min_malus (l : arbre_phylo list) : arbre_phylo = match l with
  | [] -> failwith "liste invalide d'arbres"
  | [x] -> x
  | x::suite -> min_malus_comparaison x (min_malus suite);;

(* TEST *)
let _ = assert (min_malus [arbre2; arbre] = arbre);;

(** Question 5 **)

(* 
  Fonction générant tous les arbres phylo possibles depuis 3 brins d'ADN 
  XYZ XZY
  YXZ YZX
  ZXY ZYX
*)
let gen_phylo (x : brin) (y : brin) (z : brin) : arbre_phylo list = 
  if List.length x = List.length y && List.length x = List.length z 
  then
    [br (Lf x) y (Lf z); br (Lf x) z (Lf y); br (Lf y) x (Lf z); br (Lf y) z (Lf x); br (Lf z) x (Lf y); br (Lf z) y (Lf x)]
  else 
    raise (Failure "brins de tailles différentes");;

(* TEST *)
let _ = assert (gen_phylo [A; A; T; A] [G; C; T; T] [A; A; G; A] = [
    Br (Lf [A; A; T; A], [G; C; T; T], 7, Lf [A; A; G; A]);
    Br (Lf [A; A; T; A], [A; A; G; A], 5, Lf [G; C; T; T]);
    Br (Lf [G; C; T; T], [A; A; T; A], 4, Lf [A; A; G; A]);
    Br (Lf [G; C; T; T], [A; A; G; A], 5, Lf [A; A; T; A]);
    Br (Lf [A; A; G; A], [A; A; T; A], 4, Lf [G; C; T; T]);
    Br (Lf [A; A; G; A], [G; C; T; T], 7, Lf [A; A; T; A])]);;
(* J'ai effectué un test avec un assert pour vérifier que le raise Failure marche.
    Néanmoins, le assert me soulève tout de même une Exception qui se propage.
*)

(** Question 6 **)

(*
  0 - Vérification de la taille de la liste des brins (pour construction arbre binaire complet)
      E : brin list
      S : bool

  1 - Génération de toutes les permutations possibles de la liste des brins
      E : brin list
      S : brin list list

  2 - Génération de tous les arbres possibles (1 arbre/permutation)
      E : brin list list
      S : arbre_phylo list

  3 - Détermination de l'arbre minimisant le malus global
      E : arbre_phylo list
      S : arbre_phylo
 *)

(* 
##############################################################################
                                0
##############################################################################
*)

(*
  Fonction logarithme de base b appliqué à x
 *)
let logb x b = (log10 x)/.(log10 b);;

(* Fonction de vérification si un nombre est un entier *)
let isInt x = (fst (modf x)) = 0.

(* 
  Fonction vérifiant si la taille de la liste des brins est autorisée pour
  la construction des arbres phylogénétiques
*)
let verification_taille (l : brin list) : bool = isInt (logb (float_of_int(List.length l + 1)) 2.);;


(* 
##############################################################################
                                1
##############################################################################
*)

(*
  Fonction renvoyant une liste de liste de brins où dans chaque nouvelle liste de cette liste,
  le brin elt est rajouté à chaque endroit de la liste d'origine

  Exemple : Si liste = [[A]; [T]] et elt = [C] alors
  insertion_liste [C] [[A]; [T]];; renverra [[C; A; T]; [A; C; T]; [A; T; C]]
 *)
let rec insertion_liste (elt : brin) (liste : brin list) : (brin list list) = match liste with 
  | [] -> [[elt]]
  | x::l -> (elt::liste) :: (List.map (fun e -> x::e) (insertion_liste elt l));;

(* 
  Fonction générant la liste de toutes les permutations possibles depuis une liste de brins 
  donnée en paramètre
*)
let rec generation_permutations (liste : brin list) : (brin list list) = match liste with
  | [] -> [liste]
  | x::l -> List.flatten (List.map (insertion_liste x) (generation_permutations l));;

(* 
##############################################################################
                                2
##############################################################################
*)

(* 
  Fonction permettant de découper une liste en 2 listes 
  découpées à partir du n-ième index de la liste d'origine
*)
let rec decoupe_liste (n : int) (l : brin list) =
  if n = 0 then ([], l) else
    match l with
    | [] -> ([], [])
    | h :: t -> let (l1, l2) = decoupe_liste (n-1) t in (h :: l1, l2);;

(* 
  Fonction permettant de découper une liste en 2 listes
  découpées à partir du centre de la liste d'origine
*)
let decoupe_liste_moitie (l : brin list) = decoupe_liste (List.length l / 2) l;;

(*
  Fonction retournant la première liste découpée
  (première moitiée de la liste d'origine)
*)
let get_first_half (l : brin list) : (brin list) = fst (decoupe_liste_moitie l);;

(*
  Fonction retournant la deuxième liste découpée
  (deuxième moitiée de la liste d'origine)
*)
let get_second_half (l : brin list) : (brin list) = snd (decoupe_liste_moitie l);;


(*
  Fonction de construction d'un arbre à partir d'une liste de brins
 *)
let rec generation_arbre (liste : brin list) : (arbre_phylo) = match liste with
  | [] -> Lf([])
  | [x] -> Lf(x)
  | x::l -> br (generation_arbre (get_first_half l)) x (generation_arbre (get_second_half l));;

(* TEST *)
let _ = assert (generation_arbre [[A]; [T]; [C]; [G]; [C]; [A]; [A]] = Br (Br (Lf [C], [T], 2, Lf [G]), [A], 6, Br (Lf [A], [C], 2, Lf [A])));;

(* 
  Fonction générant tous les arbres possibles depuis une liste de liste de brins
*)
let rec generation_liste_arbres (liste : brin list list) : (arbre_phylo list) = match liste with 
  | [] -> []
  | x::l -> (generation_arbre x)::(generation_liste_arbres l);;

(* 
##############################################################################
                                3
##############################################################################
*)

(* 
  Renvoie l'arbre phylogénétique qui minimise la valeur du malus global parmi tous les arbres
  que l'on peut générer
*)
let gen_min_malus_phylo (l : brin list) : arbre_phylo = 
  if verification_taille l then
    min_malus (generation_liste_arbres (generation_permutations l))
  else 
    raise (Failure "La liste ne permet pas de créer des arbres binaires parfaits.");;

(* TEST *)
let _ = assert (gen_min_malus_phylo [[A; C]; [A; T]; [A; G]] = Br (Lf [A; T], [A; G], 2, Lf [A; C]));;
let _ = assert (gen_min_malus_phylo [[A; A; T; A]; [G; C; T; A]; [A; G; C; A]; [G; A; A; C]; [C; A; G; T]; [T; A; C; A]; [G; T; G; A]] = Br (Br (Lf [G; T; G; A], [G; C; T; A], 5, Lf [G; A; A; C]), [A; A; T; A], 14,
                                                                                                                                             Br (Lf [C; A; G; T], [T; A; C; A], 5, Lf [A; G; C; A])));;

(* TEST PLUS COMPLEXES *)

(*
  Fonction
 *)
let rec arbre_taille t = match t with
  | Lf(x) -> 1
  | Br(l, x, m, r) -> 1 + arbre_taille r + arbre_taille l;;

(*
On effectue un parcours pour avoir toutes les paires de sommets en quelques sortes
ATTENTION : Cette fonction renvoie un exemplaire avec la racine en plus
 *)
let rec parcours arbre = match arbre with
  | Lf(x) -> [x]
  | Br(l, x, m, r) ->  [x] @ parcours l @ [x] @ parcours r @ [x];;

(* Il faut maintenant enlever le premier élément qui est la racine qui ne sert à rien 
   2 options :
   - la fonction suppression_racine qui enlève le premeir elt
   - la fonction List.tl
*)
(* 
  Fonction supprimant la racine d'une liste
*)
let suppression_racine (l:brin list) : (brin list) = match l with
  | [] -> []
  | x::m -> m;;

(*
  Fonction fold_pairs
 *)
let fold_pairs f v = function
  | [] 
  | [_] -> v
  | x1 :: l ->
    let acc, xn =
      List.fold_left (fun (acc, prev) x -> f acc prev x, x) (v, x1) l in
    f acc xn x1;;

(* Similitude d'un arbre *)
let sim arbre = (fold_pairs (fun acc h g -> acc +. similarite h g) 0.0 (suppression_racine (parcours arbre)));;

(* TEST *)
let _ = assert (sim (gen_min_malus_phylo [[A; A; T; A]; [G; C; T; T]; [A; A; G; A]]) = sim (min_malus (gen_phylo [A; A; T; A] [G; C; T; T] [A; A; G; A])));;