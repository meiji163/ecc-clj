# ecc-clj
A clojure implementation of Reed-Solomon and BCH codes, for educational purposes.

It also come with a QR code generator!

The associated blogpost is here: https://meiji163.github.io/post/ecc/

## Usage

``` clojure
(require '[ecc-clj.poly :as p])
(require '[ecc-clj.core :as ecc])

;; create a Galois field
(def GF256
  "GF(255) constructed as GF2[x]/<x^8+x^4+x^3+x^2+1>"
  (let [GF2-poly (p/parse-bin "100011101")]
    (p/char2-field 2 GF2-poly)))

;; create a Reed-Solomon code generating polynomial
(def RS-255-248
  "(255,248) Reed-Solomon code over GF256"
  (let [roots (take 7 (:exp GF256))]
    (reduce
     (fn [p1 p2] (p/* p1 p2 GF256))
     (for [r roots] [r 1]))))
     
;; encode a message
(def my-encoded-msg
  (let [message
        (str
         "Hello world! This is meiji163 transmitting from Neptune. "
         "It's cold here, Please send hot chocolate. Thanks. "
         "Now that I think of it, ramen would be good too if you have some.")

        data (vec (map int message))
        padded (p/shift-right
                data
                (- 223 (count data)))]
    (ecc/encode padded RS-255-223 GF256)))

;; decode with errors
(let [err [0 42 1 0 0 163 0 0 66 0 0 0 0 0 101 100]
      max-errs 8
      decode-me (p/+ my-encoded-msg err GF256)]
  (ecc/decode decode-me max-errs GF256))
;; => {:locations [1 2 5 8 14 15], 
;;     :sizes [42 1 163 66 101 100]}
```

## Test

Run tests with `lein test`

## License

Copyright © 2023 FIXME

This program and the accompanying materials are made available under the
terms of the Eclipse Public License 2.0 which is available at
http://www.eclipse.org/legal/epl-2.0.

This Source Code may also be made available under the following Secondary
Licenses when the conditions for such availability set forth in the Eclipse
Public License, v. 2.0 are satisfied: GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or (at your
option) any later version, with the GNU Classpath Exception which is available
at https://www.gnu.org/software/classpath/license.html.
