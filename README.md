# ecc-clj
A clojure implementation of Reed-Solomon and BCH codes, for educational purposes.

It also come with a QR code generator!

The associated blogpost is here: https://meiji163.github.io/post/ecc/

## Usage

``` clojure
(require '[ecc-clj.poly :as p])
(require '[ecc-clj.core :as ecc])
    
;; encode a message with (255,233) Reed-Solomon code
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
    (ecc/encode padded ecc/RS-255-223 ecc/GF256)))

;; correct message with errors
(let [err [0 42 1 0 0 163 0 0 66 0 0 0 0 0 101 100]
      max-errs 8
      decode-me (p/+ my-encoded-msg err ecc/GF256)]
  (ecc/decode decode-me max-errs ecc/GF256))
;; => {:locations [1 2 5 8 14 15], 
;;     :sizes [42 1 163 66 101 100]}


;; create a QR code
(import java.io.File)
(import java.imageio.ImageIO)
(require '[ecc-clj.qr :as qr])

(let [mask-code 1
      databits (->> "github.com/github"
                  (map int)
                  (qr/encode-bytes)
                  (qr/mask-bits mask-code))
      fmtbits (qr/fmt-bits mask-code)
      myimg (qr/qr-image fmtbits databits)]
(ImageIO/write myimg "PNG" (File. "test.png")))
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
