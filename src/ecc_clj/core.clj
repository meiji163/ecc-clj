(ns ecc-clj.core
  (:require [ecc-clj.poly :as p]))

(def GF255
  "construct GF(255) as GF2[x]/<x^8+x^7+x^2+x+1>"
  (let [GF2-poly (p/parse-bin "110000111")]
    (p/char2-field 2 GF2-poly)))

(def GF16
  "construct GF(16) as GF2[x]/<x^4+x+1>"
  (let [GF2-poly (p/parse-bin "10011")]
    (p/char2-field 2 GF2-poly)))

(def GF8
  "construct GF(8) as GF2[x]/<x^3+x+1>"
  (let [GF2-poly (p/parse-bin "1011")]
    (p/char2-field 2 GF2-poly)))

(defn encode
  "systematic encoding: the data polynomial is the
  highest k coefficients and the (n-k) lower are checks"
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

(def RS-7-5
  "(7,5) Reed-Solomon code over GF(8) with
  generating polynomial
  g(x) = (x-w)(x-w^2), where w=2 "
  (p/* [2 1] [4 1] GF8))
;; => [3 6 1]

(def RS-7-5-table
  "Decoding table for (7,5) Reed-Solomon code.
  It maps e(x) mod g(x) to e(x) for each error poly.
  The polynomials are encoded as integers"
  (let [errors (single-errors 7 8)
        rems (map #(p/mod % RS-7-5 GF8) errors)
        to-int #(p/digits-to-int % 8)]
    (zipmap
     (map to-int rems)
     (map to-int errors))))

(def RS-7-3
  "(7,3) Reed-Solomon code over GF(8)
  with generating polynomial
  g(x) = (x-w^4)(x-w^5)(x-w^6)(x-w^7), where w=2"
  (reduce
   (fn [p1 p2] (p/* p1 p2 GF8))
   [[7 1] [3 1] [6 1] [1 1]]))
;; => [7 0 5 3 1]

(def RS-7-3-table
  (let [errors (concat
                (single-errors 7 8)
                (double-errors 7 8))
        rems (map #(p/mod % RS-7-3 GF8) errors)
        to-int #(p/digits-to-int % 8)]
    (zipmap
     (map to-int rems)
     (map to-int errors))))


(defn RS-7-decode
  "decode a n=7 Reed Solomon code over GF8
  given the generating polynomial and decoding table"
  [gen-poly decoding-tbl data]
  (let [deg (dec (count gen-poly))
        extract-data (fn [v] (subvec v deg))

        rem (p/mod data gen-poly GF8)
        rem-int (p/digits-to-int rem 8)
        error (decoding-tbl rem-int)]
    (cond
      ;; no errors
      (empty? rem) (extract-data data)
      ;; uncorrectable error
      (nil? error) nil
      :else
      ;; correctable error; subtract it
      (let [err-poly (p/base-n-digits error 8)]
        (extract-data
         (p/- data err-poly GF8))))
    ))

(def RS-7-5-decode
  (partial RS-7-decode RS-7-5 RS-7-5-table))

(def RS-7-3-decode
  (partial RS-7-decode RS-7-3 RS-7-3-table))

(def RS-15-11
  "(15,11) Reed-Solomon code over GF(16)
  with generating polynomial
  g(x) = (x-w)(x-w^2)(x-w^3)(x-w^4), where w=2"
  (reduce
   (fn [p1 p2] (p/* p1 p2 GF16))
   [[2 1] [4 1] [8 1] [3 1]]))
;; => [7 8 12 13 1]

(def BCH-15-7
  "(15,7) BCH code over GF2 corrects two errors.
  g(x) = (x^4+x+1)(x^4+x^3+x^2+x+1) is constructed so that
  w,w^2,w^3,w^4 are zeroes of g(x) in GF(16)"
  (p/* [1 1 0 0 1] [1 1 1 1 1] p/binary-field))
;; => [1 0 0 0 1 0 1 1 1]

(def RS-255-223
  "(255,223) Reed-Solomon code over GF255 corrects 16 errors.
  g(x)=(x-w)(x-w^2)...(x-w^32)
  This is the standard recommended by CCSDS."
  (let [mul (:* GF255)
        prim 2
        roots (take 32
                    (iterate #(mul prim %) prim))]
    (reduce
     (fn [p1 p2] (p/* p1 p2 GF255))
     (for [r roots] [r 1]))))

(comment
  (-> [5 2 3 1 6]
      (encode RS-7-5 GF8)
      (RS-7-5-decode))
  ;; => [5 2 3 1 6]
  (encode [1 6 3] RS-7-3 GF8)
  ;; => [0 2 2 4 1 6 3]
  (RS-7-3-decode [1 2 2 4 1 7 3])
  ;; => [1 6 3]

  (encode [0 1 0 1 0 1 0] BCH-15-7 p/binary-field)
  ;; => [0 1 0 1 1 0 0 0 0 1 0 1 0 1 0]

  ;; check 2 is a primitive element
  (let [mul (:* GF255)]
    (= (range 1 256)
       (sort (take 255 (iterate #(mul 2 %) 2)))))
  )
