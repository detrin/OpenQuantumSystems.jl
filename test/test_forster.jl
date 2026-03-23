using Test
using OpenQuantumSystems

# Helper: build a simple single-mode molecule
function make_mol(E_excited; omega=200.0, hr=0.02, nvib=5)
    modes = [Mode(omega, 2*hr)]  # shift = 2*hr_factor
    Molecule(modes, nvib, [0.0, E_excited])
end

@testset "forster" begin

    @testset "absorption_spectrum" begin
        mol = make_mol(12500.0; omega=200.0, hr=0.02, nvib=5)
        ls = absorption_spectrum(mol; sigma=10.0, n_points=3000)

        # Spectrum is defined on a frequency grid
        @test length(ls.freqs) == 3000
        @test length(ls.intensities) == 3000

        # Normalized: ∫ dν ≈ 1
        dω = ls.freqs[2] - ls.freqs[1]
        @test sum(ls.intensities) * dω ≈ 1.0 atol=1e-3

        # All intensities non-negative
        @test all(ls.intensities .>= 0.0)

        # Peak should be near E_eg = 12500 (the 0-0 line dominates for small S)
        peak_freq = ls.freqs[argmax(ls.intensities)]
        @test peak_freq ≈ 12500.0 atol=50.0

        # Spectrum extends to higher frequencies (vibronic progression)
        @test ls.freqs[end] > 12500.0 + 200.0
    end

    @testset "emission_spectrum" begin
        mol = make_mol(12500.0; omega=200.0, hr=0.02, nvib=5)
        ls = emission_spectrum(mol; sigma=10.0, n_points=3000)

        dω = ls.freqs[2] - ls.freqs[1]
        @test sum(ls.intensities) * dω ≈ 1.0 atol=1e-3
        @test all(ls.intensities .>= 0.0)

        # Peak near E_eg = 12500
        peak_freq = ls.freqs[argmax(ls.intensities)]
        @test peak_freq ≈ 12500.0 atol=50.0

        # Emission extends to lower frequencies (Stokes-shifted progression)
        @test ls.freqs[1] < 12500.0 - 200.0
    end

    @testset "absorption vs emission Stokes shift" begin
        # For a displaced oscillator, absorption peaks above E_eg and
        # emission peaks below E_eg — they are mirror images around E_eg.
        mol = make_mol(12500.0; omega=200.0, hr=0.1, nvib=6)
        abs_ls = absorption_spectrum(mol; sigma=5.0, n_points=4000)
        em_ls  = emission_spectrum(mol;  sigma=5.0, n_points=4000)

        abs_mean = sum(abs_ls.freqs .* abs_ls.intensities) / sum(abs_ls.intensities)
        em_mean  = sum(em_ls.freqs  .* em_ls.intensities)  / sum(em_ls.intensities)

        # Absorption is blue-shifted relative to emission (Stokes shift)
        @test abs_mean > em_mean

        # Both means are symmetric around E_eg = 12500 (for equal modes on ground/excited)
        E_eg = 12500.0
        @test abs(abs_mean - E_eg) ≈ abs(em_mean - E_eg) atol=20.0
    end

    @testset "temperature effect on absorption" begin
        mol = make_mol(12500.0; omega=500.0, hr=0.1, nvib=5)

        abs_T0   = absorption_spectrum(mol; T=0.0,   sigma=20.0)
        abs_T300 = absorption_spectrum(mol; T=300.0, sigma=20.0)

        # Both normalized
        dω = abs_T0.freqs[2] - abs_T0.freqs[1]
        @test sum(abs_T0.intensities)   * dω ≈ 1.0 atol=1e-3
        @test sum(abs_T300.intensities) * dω ≈ 1.0 atol=1e-3

        # At T=300K, hot bands populate higher vib states, shifting mean to lower freq
        mean_T0   = sum(abs_T0.freqs   .* abs_T0.intensities)   / sum(abs_T0.intensities)
        mean_T300 = sum(abs_T300.freqs .* abs_T300.intensities) / sum(abs_T300.intensities)
        @test mean_T300 < mean_T0
    end

    @testset "spectral_overlap" begin
        mol = make_mol(12500.0; omega=200.0, hr=0.05, nvib=5)

        em_ls  = emission_spectrum(mol; sigma=30.0)
        abs_ls = absorption_spectrum(mol; sigma=30.0)

        S = spectral_overlap(em_ls, abs_ls)

        # Overlap must be positive and finite
        @test S > 0.0
        @test isfinite(S)

        # Identical donor and acceptor → non-zero overlap
        @test S > 1e-6

        # Donor and acceptor with large energy gap → small overlap
        mol_shifted = make_mol(15000.0; omega=200.0, hr=0.05, nvib=5)
        abs_shifted = absorption_spectrum(mol_shifted; sigma=30.0)
        S_mismatch = spectral_overlap(em_ls, abs_shifted)
        @test S_mismatch < S * 0.01
    end

    @testset "forster_rate" begin
        mol = make_mol(12500.0; omega=200.0, hr=0.05, nvib=5)
        J = 100.0  # cm⁻¹

        k = forster_rate(J, mol, mol)

        # Rate is positive and finite
        @test k > 0.0
        @test isfinite(k)

        # Scales as J² (double J → 4× rate)
        k2 = forster_rate(2*J, mol, mol)
        @test k2 ≈ 4 * k rtol=1e-6

        # Symmetric: k(D→A) == k(A→D) for identical molecules
        mol_a = make_mol(12500.0; omega=200.0, hr=0.05, nvib=5)
        mol_b = make_mol(12500.0; omega=200.0, hr=0.05, nvib=5)
        @test forster_rate(J, mol_a, mol_b) ≈ forster_rate(J, mol_b, mol_a)

        # Zero coupling → zero rate
        @test forster_rate(0.0, mol, mol) == 0.0

        # Far-detuned molecules → much smaller rate
        mol_far = make_mol(15000.0; omega=200.0, hr=0.05, nvib=5)
        k_far = forster_rate(J, mol, mol_far)
        @test k_far < k * 0.01
    end

    @testset "LineshapeResult fields" begin
        mol = make_mol(12000.0)
        ls = absorption_spectrum(mol; n_points=500)
        @test ls isa LineshapeResult
        @test length(ls.freqs) == 500
        @test length(ls.intensities) == 500
    end

end
