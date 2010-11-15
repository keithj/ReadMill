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

(in-package :cl-user)

(asdf:load-system :deoxybyte-systems)

(in-package :uk.co.deoxybyte-systems)

(asdf:defsystem readmill
  :name "readmill"
  :version "0.0.5"
  :author "Keith James"
  :licence "GPL v3"
  :depends-on (:deoxybyte-systems
               (:version :deoxybyte-run "0.4.6")
               (:version :cl-sam "0.10.0")
               (:version :eager-future "0.4.0")
               (:version :cl-json "0.4.0"))
  :in-order-to ((test-op (load-op :readmill :readmill-test)))
  :components ((:module :readmill
                        :serial t
                        :pathname "src/"
                        :components ((:file "package")
                                     (:file "conditions")
                                     (:file "utilities")
                                     (:file "json")
                                     (:file "filters")
                                     (:file "read-analysis")
                                     (:file "read-filtering")
                                     (:file "commands")
                                     (:file "readmill-cli")))
               (:lift-test-config :lift-tests
                                  :pathname "readmill-test"
                                  :target-system :readmill)
               (:cldoc-config :cldoc-documentation
                              :pathname "doc/html/"
                              :target-system :readmill)))
