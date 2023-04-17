package gutils

import "unicode"

func Contains[T comparable](s []T, e T) bool {
	for _, v := range s {
		if v == e {
			return true
		}
	}
	return false
}

func ContainsRune(s []uint8, e uint8, ignoreCase bool) bool {
	for _, v := range s {
		if v == e || (ignoreCase && unicode.ToLower(rune(v)) == unicode.ToLower(rune(e))) {
			return true
		}
	}
	return false
}
