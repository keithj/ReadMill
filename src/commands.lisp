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
  "Applies the base pattern report CLI to PARSED-ARGS."
  (apply #'write-pattern-report
         (mapcar (lambda (option)
                   (option-value option parsed-args))
                 '(report-file pattern-char min-freq input-file read-group))))

(defun read-filter (parsed-args &optional argv)
  "Applies the read filter CLI to PARSED-ARGS."
  (let* ((input-file (option-value 'input-file parsed-args))
         (output-file (option-value 'output-file parsed-args))
         (filter-args '(read-group min-quality queries))
         (filters (remove-if #'null
                             (mapcan (lambda (arg)
                                       (make-filters arg parsed-args))
                                     filter-args)))
         (descriptors (remove-if #'null
                                 (mapcan (lambda (arg)
                                           (make-descriptors arg parsed-args))
                                         filter-args)))
         (orphans (option-value 'orphans parsed-args))
         (json-file (option-value 'json-file parsed-args nil)))
    (filter-bam argv input-file output-file filters descriptors
               :orphans orphans  :json-file json-file)))

(defgeneric make-filters (arg parsed-args)
  (:documentation "Returns a list of filter predicates appropriate to
CLI argument ARG. The list may be empty, or contain a single element.")
  (:method (arg args)
    (declare (ignore args))
    nil)
  (:method ((arg (eql 'read-group)) args)
    (when (option-value 'read-group args nil)
      (list (complement (make-rg-p (option-value 'read-group args))))))
  (:method ((arg (eql 'min-quality)) args)
    (when (option-value 'min-quality args nil)
      (list (make-quality-p (option-value 'min-quality args)
                            :start (option-value 'read-start args 0)
                            :end (option-value 'read-end args nil)))))
  (:method ((arg (eql 'queries)) args)
    (when (option-value 'queries args nil)
      (let ((start (option-value 'start args 0))
            (end (option-value 'end args nil)))
        (mapcar (lambda (query)
                  (make-subseq-p query :start start :end end))
                (option-value 'queries args))))))

(defgeneric make-descriptors (arg parsed-args)
  (:documentation "Returns a list of filter descriptors appropriate to
CLI argument ARG. The list may be empty, or contain a single element.")
  (:method (arg args)
    (declare (ignore args))
    nil)
  (:method ((arg (eql 'read-group)) args)
    (when (option-value 'read-group args nil)
      (let ((msg (format nil "not in read-group ~s"
                         (option-value 'read-group args))))
        (list (lambda (calls true)
                (describe-filter-result calls true "read-group" msg))))))
  (:method ((arg (eql 'min-quality)) args)
    (when (option-value 'min-quality args nil)
      (let* ((start (option-value 'start args 0))
             (end (option-value 'end args nil))
             (min-quality (option-value 'min-quality args))
             (msg (format nil "a base between ~d and ~a has quality below ~d"
                          start (or end "the read end") min-quality)))
        (list (lambda (calls true)
                (describe-filter-result calls true "quality-filter" msg))))))
  (:method ((arg (eql 'queries)) args)
    (when (option-value 'queries args nil)
      (let* ((start (option-value 'start args 0))
             (end (option-value 'end args nil))
             (fmt (format nil "sequence ~~s found between ~d and ~a"
                          start (or end "the read end"))))
        (mapcar (lambda (query)
                  (lambda (calls true)
                    (describe-filter-result calls true "subseq-filter"
                                            (format nil fmt query))))
                (option-value 'queries args))))))
