#pragma once

#include <zlib.h>

#include <array>
#include <cstddef>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <streambuf>
#include <string>

namespace vcf2fna_gzip_detail {

constexpr std::size_t bufferSize = 64U * 1024U;

class inputBuffer final : public std::streambuf {
public:
	explicit inputBuffer(const char* path) : file_(gzopen(path, "rb")) {
		setg(buffer_.data(), buffer_.data(), buffer_.data());
	}

	~inputBuffer() override {
		if (file_ != nullptr) gzclose(file_);
	}

	bool isOpen() const { return file_ != nullptr; }

protected:
	int_type underflow() override {
		if (gptr() < egptr()) return traits_type::to_int_type(*gptr());
		if (file_ == nullptr) return traits_type::eof();

		const int count = gzread(file_, buffer_.data(), static_cast<unsigned int>(buffer_.size()));
		if (count < 0) {
			int errorNumber = Z_OK;
			const char* message = gzerror(file_, &errorNumber);
			throw std::runtime_error("gzip read failed: " +
				std::string(message == nullptr ? "unknown zlib error" : message));
		}
		if (count == 0) return traits_type::eof();

		setg(buffer_.data(), buffer_.data(), buffer_.data() + count);
		return traits_type::to_int_type(*gptr());
	}

private:
	gzFile file_;
	std::array<char, bufferSize> buffer_{};
};

class outputBuffer final : public std::streambuf {
public:
	outputBuffer(const char* path, std::ios_base::openmode mode) :
		file_(gzopen(path, (mode & std::ios_base::app) != 0 ? "ab" : "wb")) {
		setp(buffer_.data(), buffer_.data() + buffer_.size());
	}

	~outputBuffer() override { close(); }

	bool isOpen() const { return file_ != nullptr; }

protected:
	int_type overflow(int_type character = traits_type::eof()) override {
		if (file_ == nullptr) return traits_type::eof();
		if (traits_type::eq_int_type(character, traits_type::eof())) {
			return writeBuffer() ? traits_type::not_eof(character) : traits_type::eof();
		}
		if (pptr() == epptr() && !writeBuffer()) return traits_type::eof();
		*pptr() = traits_type::to_char_type(character);
		pbump(1);
		return traits_type::not_eof(character);
	}

	int sync() override {
		if (file_ == nullptr || !writeBuffer()) return -1;
		return gzflush(file_, Z_SYNC_FLUSH) == Z_OK ? 0 : -1;
	}

private:
	bool writeBuffer() {
		const std::ptrdiff_t pending = pptr() - pbase();
		std::ptrdiff_t written = 0;
		while (written < pending) {
			const unsigned int chunk = static_cast<unsigned int>(pending - written);
			const int result = gzwrite(file_, pbase() + written, chunk);
			if (result <= 0) return false;
			written += result;
		}
		setp(buffer_.data(), buffer_.data() + buffer_.size());
		return true;
	}

	void close() noexcept {
		if (file_ == nullptr) return;
		(void)writeBuffer();
		(void)gzclose(file_);
		file_ = nullptr;
	}

	gzFile file_;
	std::array<char, bufferSize> buffer_{};
};

} // namespace vcf2fna_gzip_detail

class igzstream final : public std::istream {
public:
	explicit igzstream(const char* path, std::ios_base::openmode = std::ios_base::in) :
		std::istream(nullptr), buffer_(path) {
		rdbuf(&buffer_);
		if (!buffer_.isOpen()) setstate(std::ios_base::failbit);
	}

private:
	vcf2fna_gzip_detail::inputBuffer buffer_;
};

class ogzstream final : public std::ostream {
public:
	explicit ogzstream(const char* path, std::ios_base::openmode mode = std::ios_base::out) :
		std::ostream(nullptr), buffer_(path, mode) {
		rdbuf(&buffer_);
		if (!buffer_.isOpen()) setstate(std::ios_base::failbit);
	}

private:
	vcf2fna_gzip_detail::outputBuffer buffer_;
};
