(ns ecc-clj.core)

(def default-field
  {:unit 1
   :zero 0
   :+ +
   :- -
   :* *
   :/ /})

(defn mod-invs
  "calculate inverses of 1..p-1 in Z/pZ
  where p is prime"
  [p]
  (let [cands (set (range 2 p))]
    (loop [n 2
           c cands
           invs {1 1}]
      (if (= n p) invs
          (let [n-inv (first
                       (filter #(= 1 (mod (* n %) p))
                               cands))]
            (recur
             (inc n)
             (disj c n-inv)
             (assoc invs n n-inv)))))
    ))

(defn prime-field [p]
  (let [inv (mod-invs p)]
    {:unit 1
     :zero 0
     :+ (fn [x y] (mod (+ x y) p))
     :- (fn [x y] (mod (- x y) p))
     :* (fn [x y] (mod (* x y) p))
     :/ (fn [x y] (mod (* x (inv y)) p))
     :inv inv}))

(def binary-field
  {:unit 1
   :zero 0
   :+ bit-xor
   :- bit-xor
   :* bit-and
   :/ bit-and})

(defn bits-to-vec
  "convert int bits into vec.
  The vector is in little-endian order"
  [b]
  (loop [n b
         bits '()]
    (if (= 0 n)
      (vec (reverse bits))
      (recur (bit-shift-right n 1)
             (conj bits (bit-and n 1))))
    ))

(defn vec-to-bits
  "convert binary vec to int"
  [v]
  (let [len (count v)]
    (loop [i 0
           acc 0]
      (if (>= i len) acc
          (recur (inc i)
                 (+ acc (bit-shift-left (v i) i))))
      )))

;; Ops for Z/2Z[x]
;; A binary polynomial is represented by its coefficients packed into bits
;; e.g. x^3 + x + 1  <--> 1011

(def bpoly-add bit-xor)
(def bpoly-sub bit-xor)

(defn bpoly-deg
  "degree of the bit-polynomial (one less than binary length)"
  [b]
  (loop [deg 0
         n 1]
    (if (> n b) (dec deg)
        (recur (inc deg) (bit-shift-left n 1)))))

(defn bpoly-mul
  [p1 p2]
  (let [deg1 (bpoly-deg p1)
        deg2 (bpoly-deg p2)]
    (loop [i 0 acc 0]
      (cond (> i deg2) acc
            (bit-test p2 i) (recur
                             (inc i)
                             (bit-xor acc (bit-shift-left p1 i)))
            :else (recur (inc i) acc)))))

(defn bpoly-mod
  "remainder of one bit-polynomial when divided by another"
  [p1 p2]
  (let [deg1 (bpoly-deg p1)
        deg2 (bpoly-deg p2)
        p2-pad (bit-shift-left p2 (- deg1 deg2))]
    (loop [i deg1
           rem p1
           div p2-pad]
      (cond (< i deg2) rem
            ;;
            (not (bit-test rem i))
            (recur (dec i)
                   rem
                   (bit-shift-right div 1))
            :else
            (recur (dec i)
                   (bit-xor rem div)
                   (bit-shift-right div 1))
            ))
    ))

(defn char2-invs
  "calculate inverses in Z/2Z[x]/<p(x)>
  given a primitive element and generating polynomial"
  [prim poly]
  (let [prim-powers
        (loop [pow 1 acc [1]]
          (let [pow-next
                (bpoly-mod (* pow prim) poly)]
            (if (= 1 pow-next) acc
                (recur
                 pow-next
                 (conj acc pow-next)))))
        non-unit (rest prim-powers)]
    (zipmap non-unit (reverse non-unit))))

(defn char2-field
  "construct field Z/2Z[x]/<p(x)>
  given a primitive element and generating polynomial"
  [prim poly]
  (let [invs (char2-invs prim poly)
        modp #(bpoly-mod % poly)]
    {:unit 1
     :zero 0
     :inv invs
     :+ bit-xor
     :- bit-xor
     :* (fn [p1 p2] (modp (bpoly-mul p1 p2)))
     :/ (fn [p1 p2]
          (modp (bpoly-mul p1 (invs p2))))}))

(defn bin-string [b]
  (Integer/toString b 2))

(defn parse-bin [s]
  (Integer/parseInt s 2))

(defn poly-mul
  "multiply two polynomials over a field"
  [p1 p2 & [field]]
  (let [field (or field default-field)
        {plus :+
         mul :*} field

        l1 (count p1)
        l2 (count p2)
        deg (dec (+ l1 l2))
        conv-j (fn [j]
                 (reduce
                  plus
                  (for [i (range (+ j 1))
                        :when (and (<= 0 (- j i))
                                   (< i l1)
                                   (< (- j i) l2))]
                    (mul (nth p1 i)
                         (nth p2 (- j i))))))]
    (vec
     (map conv-j (range deg)))))

(defn poly-scal
  "multiply polynomial by a scalar"
  [p s & [field]]
  (let [field (or field default-field)
        mul (:* field)]
    (vec (map #(mul s %) p))))

(defn monomial [n]
  (let [v (vec (repeat n 0))]
    (assoc v n 1)))

(defn poly-add [p1 p2 & [field]]
  (let [field (or field default-field)
        plus (:+ field)
        l1 (count p1)
        l2 (count p2)]
    (vec
     (for [i (range (max l1 l2))]
       (plus (nth p1 i 0) (nth p2 i 0))))))

(defn poly-sub [p1 p2 & [field]]
  (let [field (or field default-field)
        minus (:- field)
        l1 (count p1)
        l2 (count p2)]
    (vec
     (for [i (range (max l1 l2))]
       (minus (nth p1 i 0) (nth p2 i 0))))))

(defn poly-mod
  "remainder when p1 is divided by p2"
  ([p1 p2 & [field]]
   (let [field (or field default-field)
         {plus :+
          minus :-
          mul :*
          div :/
          zero :zero} field

         deg1 (dec (count p1))
         deg2 (dec (count p2))]
     (loop [rem p1]
       (let [deg (dec (count rem))
             diff (- deg deg2)]
         (cond (or (< diff 0)
                   (< deg deg2)) rem
               (= zero (last rem)) (recur (subvec rem 0 deg))
               :else
               ;; do long division step
               (let [factor
                     (poly-scal (monomial diff)
                                (div (last rem) (last p2))
                                field)
                     q (poly-mul p2 factor field)]
                 (recur (poly-sub rem q field)))))
       ))))

(defn poly-ring
  "construct polynomial ring F[x]"
  [field]
  {:unit [(:unit field)]
   :zero [(:zero field)]
   :+ (fn [p1 p2] (poly-add p1 p2 field))
   :- (fn [p1 p2] (poly-sub p1 p2 field))
   :* (fn [p1 p2] (poly-mul p1 p2 field))})

;; construct GF16 as GF2[x]/<x^4+x+1>
(def GF2-poly (parse-bin "10011"))
(def GF16 (char2-field 2 GF2-poly))

;;;
;; (15,11) Reed-Solomon code over GF(16)
;;
;; g(x) = (x-w)(x-w^2)(x-w^3)(x-w^4)
;; where w=2 is a primitive element
(def RS-poly-15-11
 (reduce
  (fn [p1 p2] (poly-mul p1 p2 GF16))
  [[2 1] [4 1] [8 1] [3 1]]))
