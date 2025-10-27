package ckalman

import (
	"github.com/jurgen-kluft/ccode/denv"
	ccore "github.com/jurgen-kluft/ccore/package"
	cunittest "github.com/jurgen-kluft/cunittest/package"
)

// GetPackage returns the package object of 'ckalman'
func GetPackage() *denv.Package {
	// Dependencies
	unittestpkg := cunittest.GetPackage()
	corepkg := ccore.GetPackage()

	// The main package
	mainpkg := denv.NewPackage("ckalman")
	mainpkg.AddPackage(unittestpkg)
	mainpkg.AddPackage(corepkg)

	// 'ckalman' library
	mainlib := denv.SetupDefaultCppLibProject("ckalman", "github.com\\jurgen-kluft\\ckalman")
	mainlib.Dependencies = append(mainlib.Dependencies, corepkg.GetMainLib())

	// 'ckalman' unittest project
	maintest := denv.SetupDefaultCppTestProject("ckalman"+"_test", "github.com\\jurgen-kluft\\ckalman")
	maintest.Dependencies = append(maintest.Dependencies, unittestpkg.GetMainLib())
	maintest.Dependencies = append(maintest.Dependencies, corepkg.GetMainLib())
	maintest.Dependencies = append(maintest.Dependencies, mainlib)

	mainpkg.AddMainLib(mainlib)
	mainpkg.AddUnittest(maintest)
	return mainpkg
}
