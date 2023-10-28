(ns ecc-clj.poly
  (:refer-clojure :exclude [+ - * mod]))

(def default-field
  {:unit 1
   :zero 0
   :+ clojure.core/+
   :- clojure.core/-
   :* clojure.core/*
   :/ /})

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
                 (clojure.core/+ acc (bit-shift-left (v i) i))))
      )))

(defn base-n-digits
  "like bits-to-vec but with arbitrary base"
  [n base]
  (loop [m n
         digs '()]
    (if (= 0 m)
      (vec (reverse digs))
      (let [q (quot m base)
            r (rem m base)]
        (recur q (conj digs r))))
    ))

(defn digits-to-int
  [digs base]
  (loop [v digs
         b 1
         acc 0]
    (if (empty? v)
      acc
      (recur (rest v)
             (* base b)
             (+ acc (* b (first v)))))
    ))

(defn bin-string [b]
  (Integer/toString b 2))

(defn parse-bin [s]
  (Integer/parseInt s 2))

(defn mod-invs
  "calculate inverses of 1..p-1 in Z/pZ where p is prime"
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

(defn scale
  "multiply polynomial by a scalar"
  [p s & [field]]
  (let [field (or field default-field)
        mul (:* field)]
    (vec (map #(mul s %) p))))

(defn monomial [n]
  (let [v (vec (repeat n 0))]
    (assoc v n 1)))

(defn + [p1 p2 & [field]]
  (let [field (or field default-field)
        plus (:+ field)
        l1 (count p1)
        l2 (count p2)]
    (vec
     (for [i (range (max l1 l2))]
       (plus (nth p1 i 0) (nth p2 i 0))))))

(defn - [p1 p2 & [field]]
  (let [field (or field default-field)
        minus (:- field)
        l1 (count p1)
        l2 (count p2)]
    (vec
     (for [i (range (max l1 l2))]
       (minus (nth p1 i 0) (nth p2 i 0))))))

(defn strip
  "strip any zeros from the right"
  [poly]
  (if (not= 0 (last poly))
    poly
    (recur (subvec poly 0
                   (dec (count poly))))))

(defn *
  "multiply two polynomials over a field"
  [p1 p2 & [field]]
  (let [field (or field default-field)
        {plus :+
         mul :*} field

        l1 (count p1)
        l2 (count p2)
        deg (dec (clojure.core/+ l1 l2))
        conv-j (fn [j]
                 (reduce
                  plus
                  (for [i (range (clojure.core/+ j 1))
                        :when (and (<= 0 (clojure.core/- j i))
                                   (< i l1)
                                   (< (clojure.core/- j i) l2))]
                    (mul (nth p1 i)
                         (nth p2 (clojure.core/- j i))))))]
    (vec
     (map conv-j (range deg)))))

(defn mod
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
             diff (clojure.core/- deg deg2)]
         (cond (or (< diff 0)
                   (< deg deg2)) (strip rem)
               (= zero (last rem)) (recur (subvec rem 0 deg))
               :else
               ;; do long division step
               (let [factor
                     (scale (monomial diff)
                            (div (last rem) (last p2))
                            field)
                     q (* p2 factor field)]
                 (recur (- rem q field)))))
       ))))

(defn ring
  "construct polynomial ring F[x]"
  [field]
  {:unit [(:unit field)]
   :zero [(:zero field)]
   :+ (fn [p1 p2] (+ p1 p2 field))
   :- (fn [p1 p2] (- p1 p2 field))
   :* (fn [p1 p2] (* p1 p2 field))})

(defn shift-right
  "shift polynomial right
  equivalent to multiplying by x^n"
  [poly n & [field]]
  (let [field (or field default-field)
        zero (:zero field)]
    (vec
     (concat (repeat n zero) poly))))


;; Ops for Z/2Z[x]
;; A binary polynomial is represented by its coefficients packed into bits
;; e.g. x^3 + x + 1  <--> 1011

(def bin+ bit-xor)
(def bin- bit-xor)

(defn bindeg
  "degree of the bit-polynomial (one less than binary length)"
  [b]
  (loop [deg 0
         n 1]
    (if (> n b) (dec deg)
        (recur (inc deg) (bit-shift-left n 1)))))

(defn bin*
  [p1 p2]
  (let [deg1 (bindeg p1)
        deg2 (bindeg p2)]
    (loop [i 0 acc 0]
      (cond (> i deg2) acc
            (bit-test p2 i) (recur
                             (inc i)
                             (bit-xor acc (bit-shift-left p1 i)))
            :else (recur (inc i) acc)))))

(defn binmod
  "remainder of one bit-polynomial when divided by another"
  [p1 p2]
  (let [deg1 (bindeg p1)
        deg2 (bindeg p2)
        p2-pad (bit-shift-left p2 (clojure.core/- deg1 deg2))]
    (loop [i deg1
           rem p1
           div p2-pad]
      (cond (< i deg2) rem
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
                (binmod (bin* pow prim) poly)]
            (if (= 1 pow-next) acc
                (recur
                 pow-next
                 (conj acc pow-next)))))
        non-unit (rest prim-powers)]
    (assoc
     (zipmap non-unit (reverse non-unit))
     ;; 1 is its own inverse
     1 1)))

(defn char2-times-table
  [poly]
  (let [n-elts (bit-shift-left 1 (bindeg poly))
        pairs (for [i (range n-elts)
                    j (range n-elts)]
                [i j])]
    (zipmap
     pairs
     (map (fn [[i j]]
            (binmod (bin* i j) poly))
          pairs))
    ))

(defn char2-field
  "construct field Z/2Z[x]/<p(x)>
  given a primitive element and generating polynomial"
  [prim poly]
  (let [invs (char2-invs prim poly)
        times-table (char2-times-table poly)]
    {:unit 1
     :zero 0
     :inv invs
     :+ bit-xor
     :- bit-xor
     :* (fn [p1 p2] (times-table [p1 p2]))
     :/ (fn [p1 p2] (times-table [p1 (invs p2)]))}
    ))
