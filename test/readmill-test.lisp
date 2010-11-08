;;
;;; Copyright (C) 2008-2010 Genome Research Ltd. All rights reserved.
;;;
;;; This file is part of readmill.
;;;
;;; This program is free software: you can redistribute it and/or modify
;;; it under the terms of the GNU General Public License as published by
;;; the Free Software Foundation, either version 3 of the License, or
;;; (at your option) any later version.
;;;
;;; This program is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;; GNU General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with this program.  If not, see <http://www.gnu.org/licenses/>.
;;;

(in-package :uk.ac.sanger.readmill-test)

(deftestsuite readmill-tests ()
  ())

(defun seq-with-n (length &rest positions)
  "Returns a sequence of LENGTH with Ns at POSITIONS and random bases
elsewhere."
  (loop
     with seq = (make-string length)
     for i from 0 below length
     for base = (if (member i positions)
                    #\N
                    (ecase (random 4)
                      (0 #\A)
                      (1 #\C)
                      (2 #\G)
                      (3 #\T)))
     do (setf (char seq i) base)
     finally (return seq)))

(addtest (readmill-tests) base-patterns/1
  (let ((bamfile (fake-bam-file
                  (dxi:tmp-pathname
                   :tmpdir (merge-pathnames "data"):type "bam")
                  :num-refs 1 :ref-length 500 :read-length 10
                  :seq-fn (lambda (length)
                            (seq-with-n length 5 7 9)))))
    (with-bam (in () bamfile)
      (let ((expected '((".....N.N.N" . 10)))
            (patterns (base-patterns in #\N #\.)))
        (ensure (equalp expected patterns)
                :report "Expected ~a but found ~a"
                :arguments (expected patterns))))
    (delete-file bamfile)))

