(ns ecc-clj.qr
  (:import java.awt.image.BufferedImage
           java.awt.Color
           javax.imageio.ImageIO
           java.io.File)
  (:require [ecc-clj.poly :as p]
            [ecc-clj.core :as ecc]))

(def GF256
  "GF(255) constructed as GF2[x]/<x^8+x^4+x^3+x^2+1>"
  (let [GF2-poly (p/parse-bin "100011101")]
    (p/char2-field 2 GF2-poly)))

;; ISO/IEC 18004 Annex A:
;; The (27,19) code is shortened from (255,248) code with this polynomial
;; w^21 + w^102 x + w^238 x^2 + w^149 x^3 + w^146 x^4 + w^229 x^5 + w^87 x^6 + x^7
(def RS-255-248
  "(255,248) Reed-Solomon code over GF256.
  g(x)=(x-w^0)(x-w^1)...(x-w^6)"
  (let [roots (take 7 (:exp GF256))]
    (reduce
     (fn [p1 p2] (p/* p1 p2 GF256))
     (for [r roots] [r 1]))))

;; https://en.wikipedia.org/wiki/QR_code#/media/File:QR_Character_Placement.svg
;;                         1                   2
;;     0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0
;;
;;  0  1 1 1 1 1 1 1 0 f         0 1 1 1 1 1 1 1
;;  1  1 0 0 0 0 0 1 0 f         0 1 0 0 0 0 0 1
;;  2  1 0 1 1 1 0 1 0 f         0 1 0 1 1 1 0 1
;;  3  1 0 1 1 1 0 1 0 f         0 1 0 1 1 1 0 1
;;  4  1 0 1 1 1 0 1 0 f         0 1 0 1 1 1 0 1
;;  5  1 0 0 0 0 0 1 0 f         0 1 0 0 0 0 0 1
;;  6  1 1 1 1 1 1 1 0 1 0 1 0 1 0 1 1 1 1 1 1 1
;;  7  0 0 0 0 0 0 0 0 f         0 0 0 0 0 0 0 0
;;  8  f f f f f f 1 f f         f f f f f f f f
;;  9              0                           x
;; 10              1                           x
;; 11              0                           x
;; 12              1                           x
;; 13  0 0 0 0 0 0 0 0 1                       x
;; 14  1 1 1 1 1 1 1 0 f                       x
;; 15  1 0 0 0 0 0 1 0 f                       x
;; 16  1 0 1 1 1 0 1 0 f                       x
;; 17  1 0 1 1 1 0 1 0 f                       x
;; 18  1 0 1 1 1 0 1 0 f                       x
;; 19  1 0 0 0 0 0 1 0 f                       x
;; 20  1 1 1 1 1 1 1 0 f                       x

(def fixed1s
  (let [sqr [[0 0] [0 1] [0 2] [0 3] [0 4] [0 5] [0 6]
             [6 0] [6 1] [6 2] [6 3] [6 4] [6 5] [6 6]
             [1 0] [2 0] [3 0] [4 0] [5 0]
             [1 6] [2 6] [3 6] [4 6] [5 6]
             ;;
             [2 2] [2 3] [2 4]
             [3 2] [3 3] [3 4]
             [4 2] [4 3] [4 4]]]
    (vec (concat
          sqr
          (for [[i j] sqr] [(+ i 14) j])
          (for [[i j] sqr] [i (+ j 14)])
          [[13 8] ;; one unused format bit
           [8 6] [10 6] [12 6]
           [6 8] [6 10] [6 12]]))))

(def fixed0s
  (let [sqr [[1 1] [1 2] [1 3] [1 4] [1 5]
             [5 1] [5 2] [5 3] [5 4] [5 5]
             [2 1] [3 1] [4 1]
             [2 5] [3 5] [4 5]
             ;;
             [0 7] [1 7] [2 7] [3 7] [4 7] [5 7] [6 7] [7 7]
             [7 0] [7 1] [7 2] [7 3] [7 4] [7 5] [7 6]]]
    (vec (concat
          sqr
          (for [[i j] sqr] [(- 20 i) j])
          (for [[i j] sqr] [i (- 20 j)])
          [[9 6] [11 6]
           [6 9] [6 11]]))))

(defn v-pattern
  "(len x 2) vertical blocks either :up or :down"
  [dir xy len]
  (let [[row col] xy
        rng (cond
              (= dir :up) (reverse
                           (range (inc (- row len)) (inc row)))
              (= dir :down) (range row (+ row len))
              :else
              (throw (IllegalArgumentException.)))]
    (->> (for [r rng] [[r col] [r (dec col)]])
         (reduce concat)
         (vec))))

(defn h-pattern
  "(2 x 8) horizonal blocks either :cwise (clockwise) or :anticwise"
  [dir xy]
  (let [[row col] xy
        r+1 (inc row) r-1 (dec row)
        c-1 (- col 1) c-2 (- col 2) c-3 (- col 3)]
    (cond (= dir :cwise)
          ;; starts at upper left corner
          [[row col] [row c-1] [r+1 col] [r+1 c-1]
           [r+1 c-2] [r+1 c-3] [row c-2] [row c-3]]

          (= dir :anticwise)
          ;; starts at lower left corner
          [[row col] [row c-1] [r-1 col] [r-1 c-1]
           [r-1 c-2] [r-1 c-3] [row c-2] [row c-3]]
          :else
          (throw (IllegalArgumentException.)))))

(def v1-bit-order
  "2d coords of data bits in the order they are filled"
  (vec (reduce concat
    [(v-pattern :up [20 20] 10) (h-pattern :anticwise [10 20])
     (v-pattern :down [11 18] 8) (h-pattern :cwise [19 18])
     (v-pattern :up [18 16] 8) (h-pattern :anticwise [10 16])
     (v-pattern :down [11 14] 8) (h-pattern :cwise [19 14])
     (v-pattern :up [18 12] 12)
     ;;
     (v-pattern :up [5 12] 4) (h-pattern :anticwise [1 12])
     (v-pattern :down [2 10] 4)
     (v-pattern :down [7 10] 14)
     ;;
     (v-pattern :up [12 8] 4) (v-pattern :down [9 5] 4)
     (v-pattern :up [12 3] 4) (v-pattern :down [9 1] 4)])))

(def v1-fmt-order
  "coords of format bits. 2 copies of (5 data bits + 10 check bits)"
  [[8 0] [8 1] [8 2] [8 3] [8 4]
   [8 5] [8 7] [8 8] [7 8] [5 8] [4 8] [3 8] [2 8] [1 8] [0 8]
   ;;
   [20 8] [19 8] [18 8] [17 8] [16 8]
   [15 8] [14 8] [8 13] [8 14] [8 15] [8 16] [8 17] [8 18] [8 19] [8 20]])

(def fmt-bitmask
  [1 0 1 0 1 0 0 0 0 0 1 0 0 1 0])

(defn fmt-bits
  [mask]
  (let [mb [0 0 1]
        ec-level [0 1]
        bits (reverse (concat ec-level mb))
        enc-bt (ecc/encode bits ecc/BCH-15-5 ecc/GF2)]
    (p/+ (reverse enc-bt) fmt-bitmask ecc/GF2)))

(def alphanum
  (let [chrs "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ $%*+-./:"]
    (zipmap chrs (range (count chrs)))))

(defn encode-bytes
  "encode 17 bytes to bit vector in the format specified by QR v1"
  [bs]
  (if (not= 17 (count bs)) (throw (IllegalArgumentException.))
      (let [tobits #(p/bits-to-vec % {:pad 8 :order :big})

            len (tobits (count bs))
            end [0 0 0 0]
            enc [0 1 0 0]
            data-bits (reduce concat (map tobits bs))
            msg-bits  (reduce concat [enc len data-bits end])

            ;; NB: have to reverse twice!
            data-poly (->> msg-bits
                           (p/bits-to-bytes)
                           (reverse))
            check-poly (p/mod (p/shift-right data-poly 7)
                              RS-255-248 GF256)
            check-bits (->> check-poly
                            (reverse)
                            (map tobits)
                            (reduce concat))]
        (concat msg-bits check-bits))))

(def mask-preds
  [(fn [i j] (even? (+ i j)))
   (fn [i j] (even? i))
   (fn [i j] (= 0 (mod (+ i j) 3)))
   (fn [i j] (= 0 (mod j 3)))
   (fn [i j] (even? (+ (quot j 3) (quot i 2))))
   (fn [i j] (= 0 (+ (mod (* i j) 2) (mod (* i j) 3))))
   (fn [i j] (even? (+ (* i j) (mod (* i j) 3))))
   (fn [i j] (even? (+ i j (mod (* i j) 3))))])

(defn mask-bits
  "apply one of the QR bitmasks on the bit sequence. The mask codes are 0-7."
  [code bits]
  (let [mask? (mask-preds code)
        coord-bits (map vector v1-bit-order bits)]
    (for [[[i j] bit] coord-bits]
      (if (mask? i j)
        (bit-xor 1 bit) bit))))

(defn fill-sqr [g size coord color]
  (let [[i j] coord]
    (.setColor g color)
    (.fill g (java.awt.Rectangle. i j size size))))

(defn qr-image
  "create 21x21 QR code in BufferedImage given the format bits"
  [fmtbits databits & [scale]]
  (if (or (not= 15 (count fmtbits))
          (not= (* 26 8) (count databits)))
    (throw (IllegalArgumentException.))

    (let [scale (or scale 10)
          size (* scale (+ 21 2))
          white Color/WHITE
          black Color/BLACK
          img (BufferedImage. size size BufferedImage/TYPE_INT_RGB)
          g (.createGraphics img)

          data-coord-bits (map vector v1-bit-order databits)
          fmt-coord-bits (map vector v1-fmt-order
                              (concat fmtbits fmtbits))
          ;; coords on java image are swapped from QR coords
          img-coord (fn [[i j]]
                      [(* scale (inc j)) (* scale (inc i))])]

      (fill-sqr g size [0 0] white)
      (doall (for [crd fixed1s]
               (fill-sqr g scale (img-coord crd) black)))
      (doall (for [crd fixed0s]
               (fill-sqr g scale (img-coord crd) white)))
      (doall (for [[crd b] data-coord-bits]
               (fill-sqr g scale (img-coord crd)
                         (if (= 0 b) white black))))
      (doall (for [[crd b] fmt-coord-bits]
               (fill-sqr g scale (img-coord crd)
                         (if (= 0 b) white black))))
      (.dispose g)
      img)))

(comment
  (let [mask-code 1
        databits (->> "github.com/github"
                      (map int)
                      (encode-bytes)
                      (mask-bits mask-code))
        fmtbits (fmt-bits mask-code)
        myimg (qr-image fmtbits databits)]
    (ImageIO/write myimg "PNG" (File. "test.png")))
  )
