(ns ecc-clj.qr-test
  (:require [ecc-clj.qr :refer :all]
            [clojure.test :refer :all]))

(deftest coord-sanity-check
  (testing "coordinates add up"
    (is
     (= 26 (/ (count v1-bit-order) 8)))
    (is
     (= 30 (count v1-fmt-order)))
    (is (= (* 21 21)
           (+ (count v1-fmt-order)
              (count fixed0s)
              (count fixed1s)
              (count v1-bit-order)))))
  (testing "coordinates don't intersect"
    (is (empty? (clojure.set/intersection
                 (set fixed0s) (set fixed1s))))
    (is (empty? (clojure.set/intersection
                 (set fixed0s) (set v1-data-order))))
    (is (empty? (clojure.set/intersection
                 (set fixed1s) (set v1-data-order))))
    ))

(deftest data-encoding-test
  (let [msg "www.wikipedia.org"
        bits (vec (encode-bytes (map int msg)))]
    (testing "length is 26 bytes"
      (is (= (* 26 8) (count bits))))
    (testing "encoding bits"
      (is (= [0 1 0 0]
             (subvec bits 0 4))))
    (testing "length bits"
      (is (= [0 0 0 1 0 0 0 1]
             (subvec bits 4 12))))
    (testing "end msg bits"
      (let [lastbyte (* 18 8)]
       (is (= [0 0 0 0]
              (subvec bits (+ 4 lastbyte) (+ 8 lastbyte))))))
    ))
