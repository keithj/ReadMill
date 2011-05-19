;;;
;;; Copyright (C) 2010 Genome Research Ltd. All rights reserved.
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

(in-package :uk.ac.sanger.readmill)

(defmacro with-underscore-translation (&body body)
  "Sets up identifier transltion between Lisp and JSON."
  `(let ((json:*lisp-identifier-name-to-json* 'lisp-to-underscore)
         (json:*json-identifier-name-to-lisp* 'underscore-to-lisp))
     ,@body))

(defun lisp-to-underscore (str)
  (string-downcase
   (if (> (length str) 2)
       (substitute-if #\_ (lambda (c)
                            (char= #\- c)) str :start 1 :end (1- (length str)))
       str)))

(defun underscore-to-lisp (str)
  (string-upcase
   (if (> (length str) 2)
       (substitute-if #\- (lambda (c)
                            (char= #\_ c)) str :start 1 :end (1- (length str)))
       str)))

(defun write-json-file (filespec object)
  (with-open-file (stream filespec :direction :output :if-exists :supersede)
    (with-underscore-translation (json:encode-json object stream))
    (terpri stream)))
