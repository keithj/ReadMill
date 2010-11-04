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

(defun about (parsed-args &optional argv)
  "Reports information about the system and exits."
  (declare (ignore argv))
  (flet ((platform-info ()
           (format *standard-output* "Common Lisp ~a version ~s on ~a~%"
                   (lisp-implementation-type) (lisp-implementation-version)
                   (machine-type)))
         (version-info ()
           (format *standard-output* "ReadMill version ~s~%~%"
                   *readmill-version*)))
    (cond ((and (option-value 'platform parsed-args)
                (not (option-value 'version parsed-args)))
           (platform-info))
          ((and (option-value 'version parsed-args)
                (not (option-value 'platform parsed-args)))
           (version-info))
          (t
           (version-info)
           (platform-info)))))

(defun quality-plot (parsed-args &optional argv)
  "Applies the mean base quality plot CLI to PARSED-ARGS."
  (declare (ignore argv))
  (apply #'write-quality-plot (mapcar (lambda (option)
                                        (option-value option parsed-args))
                                      '(plot-file input-file read-group))))

(defun pattern-report (parsed-args &optional argv)
  (declare (ignore argv))
  (apply #'write-pattern-report
         (mapcar (lambda (option)
                   (option-value option parsed-args))
                 '(report-file pattern-char min-freq input-file read-group))))

(defun quality-filter (parsed-args &optional argv)
  (apply #'quality-filter-bam argv
         (mapcar (lambda (option)
                   (option-value option parsed-args))
                 '(input-file output-file read-start read-end min-quality))))

(defun subseq-filter (parsed-args &optional argv)
  (apply #'subseq-filter-bam argv
         (mapcar (lambda (option)
                   (option-value option parsed-args))
                 '(input-file output-file read-start read-end queries))))
