(ns ecc-clj.core
  (:require [ecc-clj.poly :as p]))

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

;; construct GF(16) as GF2[x]/<x^4+x+1>
(def GF16
  (let [GF2-poly (parse-bin "10011")]
   (p/char2-field 2 GF2-poly)))

;; construct GF(8) as GF2[x]/<x^3+x^2+1>
(def GF8
  (let [GF2-poly (parse-bin "1011")]
    (p/char2-field 2 GF2-poly)))

(defn encode
  [data-poly gen-poly & [field]]
  (let [field (or field p/default-field)
        shifted (p/shift-right
                 data-poly (dec (count gen-poly)))
        rem (p/mod shifted gen-poly field)]
     (p/+ shifted rem field)))

(defn single-errors
  "all single-error error polynomials"
  [len base]
  (for [i (range len)
        d (range 1 base)]
    (p/shift-right [d] i)))

(defn double-errors
  "all double-error error polynomials"
  [len base]
  (for [i (range len)
        j (range len)
        di (range 1 base)
        dj (range 1 base)
        :when (< i j)]
    (p/strip
     (assoc (vec (repeat len 0))
            i di
            j dj))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; (7,5) Reed-Solomon code over GF(8)
;; corrects 1 error
;;
;; Generating polynomial:
;; g(x) = (x-w)(x-w^2), where w=2
(def RS-7-5
  (p/* [2 1] [4 1] GF8))
;; => [3 6 1]

(def RS-7-5-table
  "decoding table for (7,5) Reed-Solomon code.
  The polynomials are encoded as integers"
  (let [errors (single-errors 7 8)
        rems (map #(p/mod % RS-7-5 GF8) errors)
        to-int #(digits-to-int % 8)]
    (zipmap
     (map to-int rems)
     (map to-int errors))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; (7,3) Reed-Solomon code over GF(8)
;; corrects 2 errors
;;
;; Generating polynomial:
;; g(x) = (x-w^4)(x-w^5)(x-w^6)(x-w^7), where w=2
(def RS-7-3
  (reduce
   (fn [p1 p2] (p/* p1 p2 GF8))
   [[7 1] [3 1] [6 1] [1 1]]))
;; => [7 0 5 3 1]

(def RS-7-3-table
  (let [errors (concat
                (single-errors 7 8)
                (double-errors 7 8))
        rems (map #(p/mod % RS-7-3 GF8) errors)
        to-int #(digits-to-int % 8)]
    (zipmap
     (map to-int rems)
     (map to-int errors))))


(defn RS-7-decode
  "decode a n=7 Reed Solomon code
  give the generating polynomial and decoding table"
  [gen-poly decoding-tbl data]
  (let [deg (dec (count gen-poly))
        extract-data (fn [v] (subvec v deg))

        rem (p/mod data gen-poly GF8)
        rem-int (digits-to-int rem 8)
        error (decoding-tbl rem-int)]
    (cond
      (empty? rem) (extract-data data)
      (nil? error) nil
      :else
      ;; subtract the error
      (let [err-poly (base-n-digits error 8)]
        (extract-data
         (p/- data err-poly GF8))))
    ))

(def RS-7-5-decode
  (partial RS-7-decode RS-7-5 RS-7-5-table))

(def RS-7-3-decode
  (partial RS-7-decode RS-7-3 RS-7-3-table))

(comment
  (-> [5 2 3 1 6]
      (encode RS-7-5 GF8)
      (RS-7-5-decode))
  (encode [1 6 3] RS-7-3 GF8)
  ;; => [0 2 2 4 1 6 3]

  (RS-7-3-decode [1 2 2 4 1 7 3]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; (15,11) Reed-Solomon code over GF(16)
;; corrects 2 errors
;;
;; Generating polynomial:
;; g(x) = (x-w)(x-w^2)(x-w^3)(x-w^4), where w=2
(def RS-15-11
  (reduce
   (fn [p1 p2] (p/* p1 p2 GF16))
   [[2 1] [4 1] [8 1] [3 1]]))
;; => [7 8 12 13 1]
