;;;
;;; Copyright (c) 2010-2011 Genome Research Ltd. All rights reserved.
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

(defvar *readmill-name* "ReadMill")
(defvar *readmill-version* "0.1.0"
  "The version number of the ReadMill application.")

(deftype quality-score ()
  "Phred quality score."
  '(integer 0 100))

(defmacro check-readmill-arguments (test-form arguments &optional error-message
                                    &rest message-arguments)
  "Checks the validity of ARGUMENTS. If TEST-FORM returns false an
{define-condition invalid-argument-error} is raised. The default error
message may be refined with an additional ERROR-MESSAGE.

Arguments:

- test-form (form): A form to be evaluated. If the form returns NIL,
  an error is raised.
- arguments (list symbols): A list of symbols to which argument values
  are bound.

Optional:

- error-message (string): An error message string.

Rest:

- message-arguments (forms): Forms that evaluate to arguments for the
  error message."
  `(progn
     (unless ,test-form
       (error 'readmill-argument-error
              :parameters ',arguments
              :arguments (list ,@arguments)
              :format-control ,error-message
              :format-arguments (list ,@message-arguments)))
     t))

(defun maybe-standard-stream (name)
  "Returns a standard stream if NAME is STRING-EQUAL to one of
\"stdin\" or \"stdout\", otherwise returns NAME."
  (cond ((string-equal "stdin" name)
         *standard-input*)
        ((string-equal "stdout" name)
         *standard-output*)
        (t
         name)))

(defun any (&rest fns)
  "Returns a new function that takes one optional argument. When the
argument is non-null, it returns T if a call to any of FNS with that
argument would return T."
  (lambda (&optional arg)
    (and arg (some (lambda (fn)
                     (funcall fn arg)) fns))))

(defun peek-alignment (bam)
  "Returns the next alignment in BAM stream, or NIL."
  (let ((pos (bgzf-tell bam)))
    (prog1
        (read-alignment bam)
      (bgzf-seek bam pos))))

(defun peek-read-length (bam)
  "Returns the length of the next read BAM stream."
  (let ((aln (peek-alignment bam)))
    (if aln
        (read-length aln)
        (error 'invalid-operation-error
               :format-control (txt "cannot peek read length;"
                                    "reached end of stream")))))

(defun decode-phred-quality (c)
  "Returns the Phred quality score encoded by the character C."
  (- (char-code c) 33))

(defun vector-mean (sums count)
  (flet ((mean (sum)
           (float (/ sum count))))
    (let ((means (make-array (length sums) :element-type 'float
                             :initial-element 0.f0)))
      (map-into means #'mean sums))))

(defun plot-per-base (plot-filespec values sample-name y-label)
  (plot-png nil values plot-filespec :title sample-name
            :x-label "Base position" :y-label y-label))

(defun plot-png (x y png-filespec &key title x-label y-label)
  "Plots X against Y and saves the result to file PNG-FILESPEC. The
plot is annotated with TITLE, X-LABEL and Y-LABEL."
  (let ((plotter (dxr:run-gnuplot))
        (plot (make-instance
               'dxr:2d-plot
               :title title
               :x-axis (make-instance 'dxr:axis :label x-label :position :x)
               :y-axis (make-instance 'dxr:axis :label y-label :position :y)
               :series (make-instance 'dxr:xy-series
                                      :x-values (or x (iota (length y)))
                                      :y-values y
                                      :style '(:linespoints
                                               :smooth
                                               :csplines)))))
    (dxr:draw-plot plotter plot :terminal :png :output png-filespec)
    (dxr:stop-gnuplot plotter)))

(defun unambiguous-pp (header)
  "Returns the previous program to be run on the data, if such can be
deduced unambiguously. If there are more than one leaf programs in the
tree, returns NIL."
  (let ((pp (last-programs header)))
    (when (endp (rest pp))
      (first pp))))

(defun add-readmill-pg (header argv)
  "Adds a PG record to BAM list HEADER based on ReadMill invocation ARGV."
  (let ((pg (pg-record "readmill-chunk"
                       :program-name *readmill-name*
                       :program-version *readmill-version*
                       :previous-program (unambiguous-pp header)
                       :command-line (format nil "~{~a~^ ~}" argv))))
    (add-pg-record header pg)))

(defun header-string (header)
  "Returns a BAM HEADER as a string."
  (with-output-to-string (s)
    (write-sam-header header s)))
